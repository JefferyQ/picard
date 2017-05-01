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
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

public class NewIlluminaBasecallsConverter<CLUSTER_OUTPUT_RECORD> extends BasecallsConverter<CLUSTER_OUTPUT_RECORD> {
    private static final Log log = Log.getInstance(NewIlluminaBasecallsConverter.class);
    Map<Integer, File> filterFileMap = new HashMap<>();

    final private Map<String, SortingCollection<CLUSTER_OUTPUT_RECORD>> barcodeToRecordCollection =
            new HashMap<>();

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
     *                                 otherwise will throw an exception
     */
    public NewIlluminaBasecallsConverter(final File basecallsDir, File barcodesDir, final int lane,
                                         final ReadStructure readStructure,
                                         final Map<String, ? extends BasecallsConverter.ConvertedClusterDataWriter<CLUSTER_OUTPUT_RECORD>> barcodeRecordWriterMap,
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
                                         final boolean ignoreUnexpectedBarcodes) {

        super(barcodeRecordWriterMap, maxReadsInRamPerTile, tmpDirs, codecPrototype, ignoreUnexpectedBarcodes,
                demultiplex, outputRecordComparator, includeNonPfReads, bclQualityEvaluationStrategy,
                outputRecordClass, numProcessors, firstTile, tileLimit, new IlluminaDataProviderFactory(basecallsDir,
                        barcodesDir, lane, readStructure, bclQualityEvaluationStrategy));

        File laneDir = new File(basecallsDir, IlluminaFileUtil.longLaneStr(lane));

        File[] cycleDirs = IOUtil.getFilesMatchingRegexp(laneDir, IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN);

        //CBCLs
        List<File> cbcls = new ArrayList<>();
        Arrays.asList(cycleDirs).forEach(cycleDir -> {
            cbcls.addAll(Arrays.asList(IOUtil.getFilesMatchingRegexp(cycleDir, "^" + IlluminaFileUtil.longLaneStr(lane) + "_(\\d{1,5}).cbcl$")));
        });
        IOUtil.assertFilesAreReadable(cbcls);

        //locs
        File locsFile = new File(basecallsDir.getParentFile(), "s.locs");
        IOUtil.assertFileIsReadable(locsFile);
        //filter

        File[] filterFiles = getTiledFiles(laneDir, Pattern.compile(ParameterizedFileUtil.escapePeriods(
                ParameterizedFileUtil.makeLaneTileRegex(".filter", lane))));
        IOUtil.assertFilesAreReadable(Arrays.asList(filterFiles));

        this.factory.setApplyEamssFiltering(applyEamssFiltering);
    }

    private File[] getTiledFiles(File baseDirectory, Pattern pattern) {
        return IOUtil.getFilesMatchingRegexp(baseDirectory, pattern);
    }

    public void doProcessing() {
        final BaseIlluminaDataProvider dataProvider = factory.makeDataProvider(true);

        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            readProgressLogger.record(null, 0);
            // If this cluster is passing, or we do NOT want to ONLY emit passing reads, then add it to the next
            if (cluster.isPf() || includeNonPfReads) {
                final String barcode = (demultiplex ? cluster.getMatchedBarcode() : null);
                addRecord(barcode, converter.convertClusterToOutputRecord(cluster));
            }
        }

        dataProvider.close();
    }

    public synchronized void addRecord(final String barcode, final CLUSTER_OUTPUT_RECORD record) {
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
