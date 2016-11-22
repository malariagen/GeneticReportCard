package org.cggh.bam;

import org.cggh.bam.target.TargetLocus;
import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public abstract class BamConfig extends BaseConfig {
	
	protected String  propPrefix;
	protected TargetLocus[] loci;
	protected boolean remapAllReads;
	protected boolean skipUnmappedReadsAnalysis;
	
	public BamConfig (File configFile, String propPrefix) throws AnalysisException  {
		super(configFile);
		this.propPrefix = propPrefix;
		loci = parseLocusConfig ();
		skipUnmappedReadsAnalysis = true;
		for (int i = 0; i < loci.length; i++) {
			if (loci[i].getAnalyzeUnmappedReads()) {
				skipUnmappedReadsAnalysis = false;
			}
		}
		remapAllReads = this.getBooleanProperty(propPrefix+"remapAllReads", false);	
	}
	
	public abstract TargetLocus[] parseLocusConfig () throws AnalysisException;
	
	public TargetLocus[] getLoci() {
		return loci;
	}
	
	public boolean getRemapAllReads() {
		return remapAllReads;
	}
	
	public boolean getSkipUnmappedReadsAnalysis() {
		return skipUnmappedReadsAnalysis;
	}
}
