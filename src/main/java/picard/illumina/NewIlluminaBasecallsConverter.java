package picard.illumina;


import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ParameterizedFileUtil;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.LocsFileReader;
import picard.util.IlluminaUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewIlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(NewIlluminaBasecallsConverter.class);
    private final Map<String, BarcodeMetric> barcodesMetrics = new HashMap<>();
    private final BarcodeMetric noMatchMetric;
    private final List<File> cbcls;
    private final List<AbstractIlluminaPositionFileReader.PositionInfo> locs = new ArrayList<>();
    private final File[] filterFiles;
    private final int maxNoCalls;
    private final int maxMismatches;
    private final int minMismatchDelta;
    private final int minimumBaseQuality;
    private final MetricsFile<BarcodeMetric, Integer> metrics;
    private final File metricsFile;
    private final Map<String, BlockingQueue<CLUSTER_OUTPUT_RECORD>> blockingQueueMap = new HashMap<>();

    /**
     * @param basecallsDir             Where to read basecalls from.
     * @param barcodesDir              Where to read barcodes from (optional; use basecallsDir if not specified).
     * @param lane                     What lane to process.
     * @param readStructure            How to interpret each cluster.
     * @param barcodeRecordWriterMap   Map from barcode to CLUSTER_OUTPUT_RECORD writer.  If demultiplex is false, must contain
     *                                 one writer stored with key=null.
     * @param demultiplex              If true, output is split by barcode, otherwise all are written to the same output stream.
     * @param maxReadsInRamPerTile     Configures number of reads each tile will store in RAM before spilling to disk.
     * @param tmpDirs                  For SortingCollection spilling.
     * @param numProcessors            Controls number of threads.  If <= 0, the number of threads allocated is
     *                                 available cores - numProcessors.
     * @param firstTile                (For debugging) If non-null, start processing at this tile.
     * @param tileLimit                (For debugging) If non-null, process no more than this many tiles.
     * @param outputRecordComparator   For sorting output records within a single tile.
     * @param codecPrototype           For spilling output records to disk.
     * @param outputRecordClass        Inconveniently needed to create SortingCollections.
     * @param includeNonPfReads        If true, will include ALL reads (including those which do not have PF set)
     * @param ignoreUnexpectedBarcodes If true, will ignore reads whose called barcode is not found in barcodeRecordWriterMap,
     * @param maxNoCalls
     * @param maxMismatches
     * @param minMismatchDelta
     * @param minimumBaseQuality
     * @param metrics
     * @param metricsFile
     */
    public NewIlluminaBasecallsConverter(final File basecallsDir, File barcodesDir, final int lane,
                                         final ReadStructure readStructure,
                                         final Map<String, ? extends ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
                                         final boolean demultiplex,
                                         final int maxReadsInRamPerTile,
                                         final List<File> tmpDirs, final int numProcessors,
                                         final Integer firstTile,
                                         final Integer tileLimit,
                                         final Comparator<CLUSTER_OUTPUT_RECORD> outputRecordComparator,
                                         final SortingCollection.Codec<CLUSTER_OUTPUT_RECORD> codecPrototype,
                                         final Class<CLUSTER_OUTPUT_RECORD> outputRecordClass,
                                         final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                                         final boolean applyEamssFiltering, final boolean includeNonPfReads,
                                         final boolean ignoreUnexpectedBarcodes, int maxNoCalls, int maxMismatches,
                                         int minMismatchDelta, int minimumBaseQuality, MetricsFile<BarcodeMetric, Integer> metrics, File metricsFile) {

        super(barcodeRecordWriterMap, maxReadsInRamPerTile, tmpDirs, codecPrototype, ignoreUnexpectedBarcodes,
                demultiplex, outputRecordComparator, includeNonPfReads, bclQualityEvaluationStrategy,
                outputRecordClass, numProcessors, firstTile, tileLimit, new IlluminaDataProviderFactory(basecallsDir,
                        barcodesDir, lane, readStructure, bclQualityEvaluationStrategy));
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        this.tiles = new ArrayList<>();
        this.metrics = metrics;
        this.metricsFile = metricsFile;
        int numBarcodes = readStructure.sampleBarcodes.length();

        barcodeRecordWriterMap.keySet().forEach(barcode -> {
            if (barcode != null) {
                int pos = 0;
                String[] bcStrings = new String[numBarcodes];
                for (int i = 0; i < numBarcodes; i++) {
                    int endIndex = readStructure.sampleBarcodes.getDescriptorLengths()[i];
                    bcStrings[i] = barcode.substring(pos, endIndex + pos);
                    pos += endIndex;
                }
                this.barcodesMetrics.put(barcode, new BarcodeMetric(null, null, barcode, bcStrings));
                blockingQueueMap.put(barcode, new ArrayBlockingQueue<>(maxReadsInRamPerTile, true));
            } else {
                //we expect a lot more unidentified reads so make a bigger queue
                blockingQueueMap.put(null, new LinkedBlockingQueue<>());
            }

        });

        File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        //CBCLs
        cbcls = new ArrayList<>();
        Arrays.asList(cycleDirs)
                .forEach(cycleDir -> cbcls.addAll(
                        Arrays.asList(IOUtil.getFilesMatchingRegexp(
                                cycleDir, "^" + IlluminaFileUtil.longLaneStr(lane) + "_(\\d{1,5}).cbcl$"))));

        if (cbcls.size() == 0) {
            throw new PicardException("No CBCL files found.");
        }

        IOUtil.assertFilesAreReadable(cbcls);

        //locs
        File locsFile = new File(basecallsDir.getParentFile(), "s.locs");
        LocsFileReader locsFileReader = new LocsFileReader(locsFile);
        while (locsFileReader.hasNext()) {
            locs.add(locsFileReader.next());
        }
        IOUtil.assertFileIsReadable(locsFile);
        //filter

        Pattern laneTileRegex = Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane)));
        filterFiles = getTiledFiles(laneDir, laneTileRegex);
        for (File filterFile : filterFiles) {
            Matcher tileMatcher = laneTileRegex.matcher(filterFile.getName());
            if (tileMatcher.matches()) {
                tiles.add(Integer.valueOf(tileMatcher.group(1)));
            }
        }
        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));
        tiles.sort(TILE_NUMBER_COMPARATOR);

        this.factory.setApplyEamssFiltering(applyEamssFiltering);
        setTileLimits(firstTile, tileLimit);

        // Create BarcodeMetric for counting reads that don't match any barcode
        final String[] noMatchBarcode = new String[readStructure.sampleBarcodes.length()];
        int index = 0;
        for (final ReadDescriptor d : readStructure.descriptors) {
            if (d.type == ReadType.Barcode) {
                noMatchBarcode[index++] = StringUtil.repeatCharNTimes('N', d.length);
            }
        }

        this.noMatchMetric = new BarcodeMetric(null, null, IlluminaUtil.barcodeSeqsToString(noMatchBarcode), noMatchBarcode);

    }

    private File[] getTiledFiles(File baseDirectory, Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    @Override
    public void doTileProcessing() {

        //thread by surface tile
        ThreadPoolExecutor executorService = new ThreadPoolExecutorWithExceptions(numThreads);

        for (Integer tile : tiles) {
            executorService.submit(new TileProcessor(tile));
            //stagger by 30 seconds to avoid all threads hitting the same file at once. Once the thread pool is full
            //we can submit faster
            try {
                if (executorService.getMaximumPoolSize() > executorService.getActiveCount()) {
                    Thread.sleep(1000);
                }
            } catch (InterruptedException e) {
                throw new PicardException("Interrupted during submit sleep.", e);
            }
        }

        executorService.shutdown();

        ThreadPoolExecutor writerExecutor = new ThreadPoolExecutorWithExceptions(barcodeRecordWriterMap.keySet().size());
        List<RecordWriter> writers = new ArrayList<>();
        for (String barcode : barcodeRecordWriterMap.keySet()) {
            RecordWriter writer = new RecordWriter(barcode, blockingQueueMap.get(barcode));
            writers.add(writer);
            writerExecutor.submit(writer);
        }

        writerExecutor.shutdown();
        awaitThreadPoolTermination("Reading executor", executorService);
        //we are done reading.. signal to the writers that we are not adding any more records to their BlockingQueue

        for (RecordWriter writer : writers) {
            writer.signalDoneAdding();
        }

        awaitThreadPoolTermination("Writing executor", writerExecutor);
        if (metricsFile != null) {
            ExtractIlluminaBarcodes.finalizeMetrics(barcodesMetrics, noMatchMetric);

            for (final BarcodeMetric barcodeMetric : barcodesMetrics.values()) {
                metrics.addMetric(barcodeMetric);
            }
            metrics.addMetric(noMatchMetric);
            metrics.write(metricsFile);
            CloserUtil.close(metricsFile);
        }
    }

    private void awaitThreadPoolTermination(String executorName, ThreadPoolExecutor executorService) {
        try {
            while (!executorService.awaitTermination(300, TimeUnit.SECONDS)) {
                final int[] queuedReads = {0};
                blockingQueueMap.values().forEach(queue -> {
                    queuedReads[0] += queue.size();
                });
                log.info(String.format("%s waiting for job completion. Finished jobs - %d : Running jobs - %d : Queued jobs  - %d : Reads in queue - %d : Reads in unidentified queue - %d",
                        executorName, executorService.getCompletedTaskCount(), executorService.getActiveCount(),
                        executorService.getQueue().size(), queuedReads[0], blockingQueueMap.get(null).size()));
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private class RecordWriter implements Runnable {
        private final String barcode;
        private final BlockingQueue<CLUSTER_OUTPUT_RECORD> recordBlockingQueue;
        private boolean stillAdding = true;

        private RecordWriter(String barcode, BlockingQueue<CLUSTER_OUTPUT_RECORD> recordBlockingQueue) {
            this.barcode = barcode;
            this.recordBlockingQueue = recordBlockingQueue;
        }

        @Override
        public void run() {
            if (this.barcode == null) {
                //set higher priority for the undefined barcode thread since we expect the most reads
                Thread.currentThread().setPriority(Thread.currentThread().getThreadGroup().getMaxPriority());
            }

            final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
            try {
                while (stillAdding) {
                    CLUSTER_OUTPUT_RECORD rec;
                    if ((rec = recordBlockingQueue.poll()) != null) {
                        writer.write(rec);
                        writeProgressLogger.record(null, 0);
                    } else {
                        Thread.sleep(100);
                    }
                }

                //we are done adding... now drain the queue
                for (CLUSTER_OUTPUT_RECORD aRecordBlockingQueue : recordBlockingQueue) {
                    writer.write(aRecordBlockingQueue);
                    writeProgressLogger.record(null, 0);
                }

                writer.close();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        void signalDoneAdding() {
            stillAdding = false;
        }
    }

    private class TileProcessor implements Runnable {
        private final int tileNum;

        TileProcessor(int tileNum) {
            this.tileNum = tileNum;
        }

        @Override
        public void run() {
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(cbcls, locs, filterFiles, tileNum, barcodesMetrics,
                    noMatchMetric, maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
            }

            dataProvider.close();

            log.info("Finished processing tile " + tileNum);
        }

        private void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            try {
                blockingQueueMap.get(barcode).put(record);
            } catch (InterruptedException e) {
                throw new PicardException(e.getMessage(), e);
            }
        }
    }


    private class ThreadPoolExecutorWithExceptions extends ThreadPoolExecutor {
        ThreadPoolExecutorWithExceptions(int threads) {
            super(threads, threads, 0, TimeUnit.SECONDS, new LinkedBlockingDeque<>());
        }

        @Override
        protected void afterExecute(Runnable r, Throwable t) {
            if (t == null && r instanceof Future<?>) {
                try {
                    Future<?> future = (Future<?>) r;
                    if (future.isDone()) {
                        future.get();
                    }
                } catch (CancellationException ce) {
                    t = ce;
                } catch (ExecutionException ee) {
                    t = ee.getCause();
                } catch (InterruptedException ie) {
                    Thread.currentThread().interrupt(); // ignore/reset
                }
            }
            if (t != null) {
                throw new PicardException(t.getMessage(), t);
            }
        }
    }
}
