package org.cggh.bam.codon;

import org.cggh.bam.target.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public class CodonConfig extends AlignmentTargetAnalysisConfig {
	
	public static final String PROP_PREFIX = "codon.";
	
	public CodonConfig (File configFile) throws AnalysisException  {
		super(configFile, PROP_PREFIX, true);
	}
}
