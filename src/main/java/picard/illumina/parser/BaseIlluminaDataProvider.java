package picard.illumina.parser;

import java.util.Iterator;

public interface BaseIlluminaDataProvider extends Iterator<ClusterData>, Iterable<ClusterData> {
    void close();

    void seekToTile(int seekAfterFirstRead);
}
