package org.cggh.bam.sampleClass;

import org.cggh.common.exceptions.*;
import java.io.*;

public class SampleClassConfig extends AlleleClassAnalysisConfig {
	
	public static final String PROP_PREFIX = "sampleClass.";
	
	protected String[]           classes;
	
	public SampleClassConfig (File configFile) throws AnalysisException  {
		super(configFile, PROP_PREFIX, false);
		classes = getStringListProperty(propPrefix+"classes");
	}
	
	public String[] getClasses() {
		return classes;
	}
}
