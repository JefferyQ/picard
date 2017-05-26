package picard.sam;

import picard.cmdline.CommandLineProgramTest;
import picard.sam.ValidateSamFile;

public class ValidateSamFileTest extends CommandLineProgramTest {

    @Override
    public String getCommandLineProgramName() {
        return ValidateSamFile.class.getSimpleName();
    }
}
