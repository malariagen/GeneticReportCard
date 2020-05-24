package org.cggh.bam;

import org.cggh.common.exceptions.*;
//import org.cggh.common.genome.*;
import java.io.*;


public abstract class SampleAnalysis extends BaseAnalysis {
	
	public SampleAnalysis (File refFastaFile, File outRootFolder) throws AnalysisException  {
		super(outRootFolder);
		
		if (refFastaFile != null) {
			// Load up the reference genome sequences
			ReferenceGenome.initialize(refFastaFile);
		}
	}

	public SampleAnalysis (File outRootFolder) throws AnalysisException  {
		this(null, outRootFolder);
	}


	public abstract void analyzeSample(Sample sample) throws AnalysisException;
}
