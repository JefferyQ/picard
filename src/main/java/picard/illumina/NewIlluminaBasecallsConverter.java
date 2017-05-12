package picard.illumina;


import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.illumina.parser.BaseIlluminaDataProvider;
import picard.illumina.parser.ClusterData;
import picard.illumina.parser.IlluminaDataProviderFactory;
import picard.illumina.parser.IlluminaFileUtil;
import picard.illumina.parser.ParameterizedFileUtil;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.LocsFileReader;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewIlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(NewIlluminaBasecallsConverter.class);
    private static final long FIFTEEN_SECONDS = 15 * 1000;
    private final Map<String, BarcodeMetric> barcodesMetrics = new HashMap<>();
    private final List<File> cbcls;
    private final List<AbstractIlluminaPositionFileReader.PositionInfo> locs = new ArrayList<>();
    private final File[] filterFiles;
    private final int maxNoCalls;
    private final int maxMismatches;
    private final int minMismatchDelta;
    private final int minimumBaseQuality;
    private final Map<String, List<SortingCollection<CLUSTER_OUTPUT_RECORD>>> barcodeRecords = new HashMap<>();

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
                                         final boolean ignoreUnexpectedBarcodes, int maxNoCalls, int maxMismatches, int minMismatchDelta, int minimumBaseQuality) {

        super(barcodeRecordWriterMap, maxReadsInRamPerTile, tmpDirs, codecPrototype, ignoreUnexpectedBarcodes,
                demultiplex, outputRecordComparator, includeNonPfReads, bclQualityEvaluationStrategy,
                outputRecordClass, numProcessors, firstTile, tileLimit, new IlluminaDataProviderFactory(basecallsDir,
                        barcodesDir, lane, readStructure, bclQualityEvaluationStrategy));
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        this.tiles = new ArrayList<>();

        barcodeRecordWriterMap.keySet().forEach(barcode -> {
            if (barcode != null)
                this.barcodesMetrics.put(barcode, new BarcodeMetric(barcode, barcode, barcode, new String[]{barcode}));
        });

        File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        //CBCLs
        cbcls = new ArrayList<>();
        Arrays.asList(cycleDirs).forEach(cycleDir -> {
            cbcls.addAll(Arrays.asList(IOUtil.getFilesMatchingRegexp(cycleDir, "^" + IlluminaFileUtil.longLaneStr(lane) + "_(\\d{1,5}).cbcl$")));
        });
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

        this.factory.setApplyEamssFiltering(applyEamssFiltering);
    }

    private File[] getTiledFiles(File baseDirectory, Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    void doTileProcessing() {

        //thread by surface tile
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.SECONDS, new LinkedBlockingDeque<>()) {
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
        };

        for (Integer tile : tiles) {
            executorService.submit(new TileProcessor(tile));
            //stagger by 30 seconds to avoid all threads hitting the same file at once.
            try {
                Thread.sleep(FIFTEEN_SECONDS);
            } catch (InterruptedException e) {
                throw new PicardException("Interrupted during submit sleep.", e);
            }
        }

        executorService.shutdown();

        try {
            while (!executorService.awaitTermination(300, TimeUnit.SECONDS)) {
                log.info(String.format("Waiting for job completion. Finished jobs - %d : Running jobs - %d : Queued jobs  - %d",
                        executorService.getCompletedTaskCount(), executorService.getActiveCount(), executorService.getQueue().size()));
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        ExecutorService writerExecutor = Executors.newFixedThreadPool(barcodeRecordWriterMap.keySet().size());

        for (String barcode : barcodeRecordWriterMap.keySet()) {
            writerExecutor.submit(new RecordWriter(barcode, this.barcodeRecords.get(barcode)));
        }
    }

    private class RecordWriter implements Runnable {
        private final String barcode;
        private final List<SortingCollection<CLUSTER_OUTPUT_RECORD>> recordsList;

        RecordWriter(String barcode, List<SortingCollection<CLUSTER_OUTPUT_RECORD>> recordsList) {
            this.barcode = barcode;
            this.recordsList = recordsList;
        }

        @Override
        public void run() {
            final ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD> writer = barcodeRecordWriterMap.get(barcode);
            for (SortingCollection<CLUSTER_OUTPUT_RECORD> records : recordsList) {
                for (CLUSTER_OUTPUT_RECORD rec : records) {
                    writer.write(rec);
                    writeProgressLogger.record(null, 0);
                }
            }
            log.debug(String.format("Closing file for barcode %s.", barcode));
            writer.close();
        }
    }


    private class TileProcessor implements Runnable {
        private final int tileNum;
        final private Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection =
                new HashMap<>();

        TileProcessor(int tileNum) {
            this.tileNum = tileNum;
        }

        @Override
        public void run() {
            final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(cbcls, locs, filterFiles, tileNum, barcodesMetrics,
                    maxNoCalls, maxMismatches, minMismatchDelta, minimumBaseQuality);

            while (dataProvider.hasNext()) {
                final ClusterData cluster = dataProvider.next();
                readProgressLogger.record(null, 0);
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
            }

            dataProvider.close();
            //we are done adding records
            this.barcodeToRecordCollection.values().forEach(SortingCollection::doneAdding);

            for (Map.Entry<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> entry : barcodeToRecordCollection.entrySet()) {
                if (barcodeRecords.containsKey(entry.getKey())) {
                    barcodeRecords.get(entry.getKey()).add(entry.getValue());
                } else {
                    List<SortingCollection<CLUSTER_OUTPUT_RECORD>> collectionList = new ArrayList<>();
                    collectionList.add(entry.getValue());
                    barcodeRecords.put(entry.getKey(), collectionList);
                }
            }
            log.info("Finished processing tile " + tileNum);
        }

        private synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
            // Grab the existing collection, or initialize it if it doesn't yet exist
            SortingCollection<CLUSTER_OUTPUT_RECORD> recordCollection = this.barcodeToRecordCollection.get(barcode);
            if (recordCollection == null) {
                // TODO: The implementation here for supporting ignoreUnexpectedBarcodes is not efficient,
                // but the alternative is an extensive rewrite.  We are living with the inefficiency for
                // this special case for the time being.
                if (!barcodeRecordWriterMap.containsKey(barcode)) {
                    if (ignoreUnexpectedBarcodes) {
                        return;
                    }
                    throw new PicardException(String.format("Read records with barcode %s, but this barcode was not expected.  (Is it referenced in the parameters file?)", barcode));
                }
                recordCollection = newSortingCollection();
                this.barcodeToRecordCollection.put(barcode, recordCollection);
            }
            recordCollection.add(record);
        }

        private synchronized SortingCollection<CLUSTER_OUTPUT_RECORD> newSortingCollection() {
            final int maxRecordsInRam =
                    Math.max(1, maxReadsInRamPerTile /
                            barcodeRecordWriterMap.size());
            return SortingCollection.newInstance(
                    outputRecordClass,
                    codecPrototype.clone(),
                    outputRecordComparator,
                    maxRecordsInRam,
                    tmpDirs);
        }
    }


}
