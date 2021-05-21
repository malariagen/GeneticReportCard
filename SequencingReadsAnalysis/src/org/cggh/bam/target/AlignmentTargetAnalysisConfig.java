package org.cggh.bam.target;

import org.cggh.common.exceptions.AnalysisException;
import java.io.File;


public class AlignmentTargetAnalysisConfig extends TargetAnalysisConfig {

	public AlignmentTargetAnalysisConfig(File configFile, String propPrefix, boolean useBamAlignment) throws AnalysisException {
		super(configFile, propPrefix, useBamAlignment);
	}

	@Override
	protected Target createTarget (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		return new AlignmentTarget (name, targetCoords, isReverse);
	}
	
}