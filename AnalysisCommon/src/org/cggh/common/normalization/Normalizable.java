package org.cggh.common.normalization;

public interface Normalizable {
	
	public boolean useRawValueForDistribution();
	
	public double getRawValue();
	
	public void setNormalizedValue(double normalizedValue);
	
}

