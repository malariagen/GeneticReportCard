package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.target.TargetLocus;
import org.cggh.bam.target.alleleClasses.AlleleClassLocus;
import org.cggh.common.exceptions.*;
import java.io.*;

public class SampleClassConfig extends BamConfig {
	
	public static final String PROP_PREFIX = "sampleClass.";
	
	protected String[]           classes;
	
	
	public SampleClassConfig (File configFile) throws AnalysisException  {
		super(configFile, PROP_PREFIX);
		classes = getStringListProperty(propPrefix+"classes");
	}
	
	public AlleleClassLocus[] getLoci() {
		return (AlleleClassLocus[])loci;
	}
	
	public String[] getClasses() {
		return classes;
	}
	
	@Override
	public TargetLocus[] parseLocusConfig() throws AnalysisException {
		return AlleleClassLocus.parseLocusConfig(this.configProperties, propPrefix);
	}
}
