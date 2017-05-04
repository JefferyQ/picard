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
    private final Iterator<AbstractIlluminaPositionFileReader.PositionInfo> locs;
    private final CbclReader reader;
    private final ReadStructure outputReadStructure;
    private final boolean usingQualityScores;
    private final Map<String, BarcodeMetric> barcodeMetricMap;

    public NewIlluminaDataProvider(final List<File> cbcls, final List<AbstractIlluminaPositionFileReader.PositionInfo> locs,
                                   final File[] filterFiles, final int lane, final OutputMapping outputMapping,
                                   boolean usingQualityScores, final Map<String, BarcodeMetric> barcodesMetrics) {
        super(lane, outputMapping);
        this.locs = locs.iterator();
        this.usingQualityScores = usingQualityScores;
        this.barcodeMetricMap = barcodesMetrics;
        Map<Integer, File> filterFileMap = new HashMap<>();
        for (File filterFile : filterFiles) {
            filterFileMap.put(fileToTile(filterFile.getName()), filterFile);
        }
        this.reader = new CbclReader(cbcls, filterFileMap, outputMapping.getOutputReadLengths());
        this.outputReadStructure = outputMapping.getOutputReadStructure();
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
        AbstractIlluminaPositionFileReader.PositionInfo positionInfo = locs.next();

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
        final boolean passingFilter = cluster.isPf();

        //do we want metrics?
        final BarcodeMatch match = BarcodeMatch.findBestBarcodeAndUpdateMetrics(barcodeSubsequences, qualityScores,
                passingFilter, barcodeMetricMap, new BarcodeMetric(), 1, 1, 1, 0);
        if (match.isMatched()) {
            cluster.setMatchedBarcode(match.getBarcode());
        }
        return cluster;
    }

    protected Integer fileToTile(final String fileName) {
        final Matcher matcher = Pattern.compile("^s_(\\d+)_(\\d{1,5}).+").matcher(fileName);
        if (!matcher.matches()) {
            return null;
        }
        return Integer.parseInt(matcher.group(2));
    }
}
