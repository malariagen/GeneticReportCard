package org.cggh.common.normalization;

import java.util.ArrayList;

import org.cggh.common.util.Statistics;


public class Normalizer {

	private double meanY;
	private double stdevY;
	private ArrayList<Normalizable> values;
	
	public Normalizer () {
		this.values = new ArrayList<Normalizable>();
	}
	
	public void addValue(Normalizable value) {
		if (value != null) {
			values.add(value);
		}
	}
	
	public void normalize() {
		
		if (values.isEmpty()) {
			return;
		}
		
		// Get the raw values
		int allValuesCount = values.size();
		int distributionValuesCount = 0;
		for (int i = 0; i < allValuesCount; i++) {
			if (values.get(i).useRawValueForDistribution()) {
				distributionValuesCount++;
			}
		}

		double[] y = new double[distributionValuesCount];
		int idx = 0;
		for (int i = 0; i < allValuesCount; i++) {
			Normalizable result = values.get(i);
			if (result.useRawValueForDistribution()) {
				y[idx++] = result.getRawValue();				
			}
		}
		
		// Work out the distribution parameters
		Statistics stats = new Statistics(y);
		meanY = stats.getMean();
		stdevY = stats.getStdev();
		
		// Normalize
		for (int i = 0; i < allValuesCount; i++) {
			Normalizable value = values.get(i);
			double rawY = value.getRawValue();
			double normalizedY = (stdevY > 0.0) ? ((rawY - meanY) / stdevY) : 0.0;
			value.setNormalizedValue(normalizedY);
		}
	}

	public int getValueCount() {
		return values.size();
	}

	public double getMeanY() {
		return meanY;
	}

	public double getStdevY() {
		return stdevY;
	}
}
