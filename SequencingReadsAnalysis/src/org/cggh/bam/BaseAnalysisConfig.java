package org.cggh.bam;

import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public class BaseAnalysisConfig extends BaseConfig {
	
	public static final int    DEFAULT_MIN_CALL_READS    = 5;
	public static final int    DEFAULT_MIN_ALLELE_READS  = 2;
	public static final double DEFAULT_MIN_ALLELE_PROP   = 0.1;
	public static final int    DEFAULT_MIN_BASE_Q_SCORE  = 10;
	
    public static final String PROP_MIN_CALL_READS      = "genotype.minCallReadCount";
	public static final String PROP_MIN_ALLELE_READS    = "genotype.minAlleleReadCount";
	public static final String PROP_MIN_ALLELE_PROP     = "genotype.minAlleleReadProp";
	public static final String PROP_MIN_BASE_Q_SCORE    = "genotype.minBaseQScore";
	
	protected String  propPrefix;

	protected int     minCallReadCount;
	protected int     minAlleleReadCount;
	protected double  minAlleleReadProp;
	protected int     minBaseQScore;
	
	public BaseAnalysisConfig (File configFile, String propPrefix) throws AnalysisException  {
		super(configFile);
		this.propPrefix = propPrefix;
		
		minCallReadCount   = this.getIntProperty(propPrefix+PROP_MIN_CALL_READS,     DEFAULT_MIN_CALL_READS);	
		minAlleleReadCount = this.getIntProperty(propPrefix+PROP_MIN_ALLELE_READS,   DEFAULT_MIN_ALLELE_READS);	
		minAlleleReadProp  = this.getDoubleProperty(propPrefix+PROP_MIN_ALLELE_PROP, DEFAULT_MIN_ALLELE_PROP);	
		minBaseQScore      = this.getIntProperty(propPrefix+PROP_MIN_BASE_Q_SCORE,   DEFAULT_MIN_BASE_Q_SCORE);
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
	
	public int getMinBaseQScore() {
		return minBaseQScore;
	}
	
	public String getPrintableDisplay() {
		return "minCallReadCount = "  + getMinCallReadCount() +
		     "\nminAlleleReadCount = "+ getMinAlleleReadCount() +
             "\nminAlleleReadProp = " + getMinAlleleReadProp() +
             "\nminBaseQScore = " + getMinBaseQScore();    		
	}

}
