package org.cggh.common.subsampling;

import org.cggh.common.config.BaseConfig;
import org.cggh.common.exceptions.AnalysisException;

import java.util.Properties;

public class SubsamplingConfig extends BaseConfig {

	public static final String SUBSAMPLING_ITERATIONS_PROP   = ".subsampling.iterations";
	public static final String SUBSAMPLING_SAMPLE_COUNT_PROP = ".subsampling.sampleCount";
	
	public static final int    SUBSAMPLING_ITERATIONS_DEFAULT = 1000;
	public static final String SUBSAMPLING_SAMPLE_COUNT_AUTO = "auto";
	
	
	private String prefix;
	private int    subsamplingIterationCount;
	private int    subsamplingSampleCount;

	public SubsamplingConfig (Properties configProperties, String prefix) throws AnalysisException {
		super(configProperties);
		this.prefix = prefix;
		
		subsamplingIterationCount = getIntProperty (prefix+SUBSAMPLING_ITERATIONS_PROP, SUBSAMPLING_ITERATIONS_DEFAULT);
		
		String subCountStr = getProperty (prefix+SUBSAMPLING_SAMPLE_COUNT_PROP);
		if (SUBSAMPLING_SAMPLE_COUNT_AUTO.equals(subCountStr)) {
			subsamplingSampleCount = 0;
		} else {
			subsamplingSampleCount = getIntProperty (SUBSAMPLING_SAMPLE_COUNT_PROP, 0);
		}
	}

	public int getSubsamplingIterationCount() {
		return subsamplingIterationCount;
	}

	public int getSubsamplingSampleCount() {
		return subsamplingSampleCount;
	}
	
	public boolean isSubsamplingAuto() {
		return (subsamplingSampleCount == 0);
	}
	
		
	public String getPrintableDisplay() {
		return (
		"-- Subsampling Configuration ("+prefix+") --" +
		"\n  Subsampling iterations: " + getSubsamplingIterationCount() +
        "\n  Subsampling sample size: " + (isSubsamplingAuto() ? "auto" : getSubsamplingSampleCount()));
	}
}
