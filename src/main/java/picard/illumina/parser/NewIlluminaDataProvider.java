package picard.illumina.parser;


import picard.illumina.BarcodeMatch;
import picard.illumina.BarcodeMetric;
import picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import picard.illumina.parser.readers.CbclReader;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewIlluminaDataProvider extends BaseIlluminaDataProvider {
    private final CbclReader reader;
    private final ReadStructure outputReadStructure;
    private final boolean usingQualityScores;
    private final Map<String, BarcodeMetric> barcodeMetricMap;
    private final int maxNoCalls;
    private final int maxMismatches;
    private final int minMismatchDelta;
    private final int minimumBaseQuality;
    private final Iterator<AbstractIlluminaPositionFileReader.PositionInfo> locsIterator;

    NewIlluminaDataProvider(final List<File> cbcls, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs,
                            final File[] filterFiles, final int lane, final int tileNum, final OutputMapping outputMapping,
                            boolean usingQualityScores, final Map<String, BarcodeMetric> barcodesMetrics,
                            int maxNoCalls, int maxMismatches, int minMismatchDelta, int minimumBaseQuality) {
        super(lane, outputMapping);
        this.usingQualityScores = usingQualityScores;
        this.barcodeMetricMap = barcodesMetrics;
        this.maxNoCalls = maxNoCalls;
        this.maxMismatches = maxMismatches;
        this.minMismatchDelta = minMismatchDelta;
        this.minimumBaseQuality = minimumBaseQuality;
        Map<Integer, File> filterFileMap = new HashMap<>();
        for (File filterFile : filterFiles) {
            filterFileMap.put(fileToTile(filterFile.getName()), filterFile);
        }
        this.reader = new CbclReader(cbcls, filterFileMap, outputMapping.getOutputReadLengths(), tileNum);
        this.outputReadStructure = outputMapping.getOutputReadStructure();
        this.locsIterator = locs.iterator();
    }

    @Override
    public void close() {
        reader.close();
    }

    @Override
    public void seekToTile(int seekAfterFirstRead) {

    }

    @Override
    public boolean hasNext() {
        return reader.hasNext();
    }

    @Override
    public ClusterData next() {
        CbclData cbclData = reader.next();

        if (cbclData == null) return null;

        AbstractIlluminaPositionFileReader.PositionInfo positionInfo = locsIterator.next();

        final ClusterData cluster = new ClusterData(outputReadTypes);
        cluster.setLane(lane);
        cluster.setTile(cbclData.getTile());

        final int[] barcodeIndices = outputReadStructure.sampleBarcodes.getIndices();
        addReadData(cluster, numReads, cbclData, positionInfo);

        final byte[][] barcodeSubsequences = new byte[barcodeIndices.length][];
        final byte[][] qualityScores = usingQualityScores ? new byte[barcodeIndices.length][] : null;
        for (int i = 0; i < barcodeIndices.length; i++) {
            barcodeSubsequences[i] = cluster.getRead(barcodeIndices[i]).getBases();
            if (usingQualityScores) qualityScores[i] = cluster.getRead(barcodeIndices[i]).getQualities();
        }

        //do we want metrics?
        final BarcodeMatch match = BarcodeMatch.findBestBarcodeAndUpdateMetrics(barcodeSubsequences, qualityScores,
                true, barcodeMetricMap, new BarcodeMetric(), maxNoCalls, maxMismatches,
                minMismatchDelta, minimumBaseQuality);

        if (match.isMatched()) {
            cluster.setMatchedBarcode(match.getBarcode());
        }
        return cluster;
    }

    private Integer fileToTile(final String fileName) {
        final Matcher matcher = Pattern.compile("^s_\\d+_(\\d{1,5}).+").matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return Integer.parseInt(matcher.group(1));
    }
}
