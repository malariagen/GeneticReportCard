package org.cggh.bam;

import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public abstract class BaseAnalysisConfig extends BaseConfig {
	
	public static final int DEFAULT_MAX_READ_MISMATCHES = 10;
	
    public static final String PROP_MIN_CALL_READS      = "genotype.minCallReadCount";
	public static final String PROP_MIN_ALLELE_READS    = "genotype.minAlleleReadCount";
	public static final String PROP_MIN_ALLELE_PROP     = "genotype.minAlleleReadProp";
	public static final String PROP_MAX_INDEL_SIZE      = "alignment.maxIndelSize";
	public static final String PROP_MAX_READ_MISMATCHES = "alignment.maxReadMismatches";
	
	protected String  propPrefix;
	protected Locus[] loci;

	protected int     minCallReadCount;
	protected int     minAlleleReadCount;
	protected double  minAlleleReadProp;
	
	protected int     maxIndelSize;
	protected int     maxReadMismatches;
	protected boolean useBamAlignment;
	protected boolean analyzeUnmappedReads;
	
	public BaseAnalysisConfig (File configFile, String propPrefix, boolean useBamAlignment) throws AnalysisException  {
		super(configFile);
		this.propPrefix = propPrefix;
		
		minCallReadCount   = this.getIntProperty(propPrefix+PROP_MIN_CALL_READS,      5);	
		minAlleleReadCount = this.getIntProperty(propPrefix+PROP_MIN_ALLELE_READS,    2);	
		minAlleleReadProp  = this.getDoubleProperty(propPrefix+PROP_MIN_ALLELE_PROP,  0.1);
		maxIndelSize       = this.getIntProperty(propPrefix+PROP_MAX_INDEL_SIZE,      0);	
		maxReadMismatches  = this.getIntProperty(propPrefix+PROP_MAX_READ_MISMATCHES, DEFAULT_MAX_READ_MISMATCHES);	
		
		this.useBamAlignment = useBamAlignment;
		loci = parseLocusConfig ();
		analyzeUnmappedReads = false;
		for (int i = 0; i < loci.length; i++) {
			if (loci[i].getAnalyzeUnmappedReads()) {
				analyzeUnmappedReads = true;
				break;
			}
		}
	}
	
	public int getMinCallReadCount() {
		return minCallReadCount;
	}
	
	public int getMinAlleleReadCount() {
		return minAlleleReadCount;
	}

	public double getMinAlleleReadProp() {
		return minAlleleReadProp;
	}

	public int getMaxIndelSize() {
		return maxIndelSize;
	}

	public int getMaxReadMismatches() {
		return maxReadMismatches;
	}
	
	public boolean getUseBamAlignment () {
		return useBamAlignment;
	}
	
	public Locus[] getLoci() {
		return loci;
	}
	
	public boolean getAnalyzeUnmappedReads() {
		return analyzeUnmappedReads;
	}
	
	public abstract Locus[] parseLocusConfig () throws AnalysisException;
	
}
