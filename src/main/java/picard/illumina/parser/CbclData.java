package picard.illumina.parser;

public class CbclData extends BclData implements PfData {
    private final int tile;

    public CbclData(int[] outputLengths, int tile) {
        super(outputLengths);
        this.tile = tile;
    }

    //CBCLs currently only contain PF reads.
    @Override
    public boolean isPf() {
        return true;
    }

    public int getTile() {
        return tile;
    }
}
