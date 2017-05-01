package picard.illumina.parser;


import java.util.Iterator;

public class NewIlluminaDataProvider implements BaseIlluminaDataProvider {
    @Override
    public void close() {

    }

    @Override
    public void seekToTile(int seekAfterFirstRead) {

    }

    @Override
    public Iterator<ClusterData> iterator() {
        return null;
    }

    @Override
    public boolean hasNext() {
        return false;
    }

    @Override
    public ClusterData next() {
        return null;
    }
}
