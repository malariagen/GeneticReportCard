package org.cggh.common.util;

import java.util.ArrayList;
import java.util.Arrays;

public class Statistics {
	
	private double   min;
	private double   max;
	private double   mean;
	private double   sum;
	private double   stdev;
	private int      count;
	private double   median;
	private double   coefVar;
	private double   quartile25;
	private double   quartile75;
	private double   quantile10;
	private double   quantile90;
	private double   quantile05;
	private double   quantile95;
	
	public Statistics (double[] inputValues) {
		double[] set = Arrays.copyOf(inputValues, inputValues.length);
		processValues(set);
	}
	
	public Statistics (double[] inputValues, int fromIdx, int toIdx) {
		double[] set = Arrays.copyOfRange(inputValues, fromIdx, toIdx);
		processValues(set);
	}
	
	public Statistics (ArrayList<Double> inputValuesList) {
		int len = inputValuesList.size();		
		double[] set = new double[len];
		for (int i = 0; i < len; i++) {
			set[i] = (double)inputValuesList.get(i);
		}
		processValues(set);
	}
	
	public Statistics (int[] inputValues) {
		int len = inputValues.length;
		double[] set = new double[len];
		for (int i = 0; i < len; i++) {
			set[i] = (double)inputValues[i];
		}
		processValues(set);
	}
	
	private void invalidate() {
		count = 0;
		min = max = Double.NaN;
		mean = median = stdev = coefVar = Double.NaN;
		sum = Double.NaN;
		min = Double.NaN;
		min = Double.NaN;
		quartile25 = quartile75 = Double.NaN;
		quantile10 = quantile90 = Double.NaN;
		quantile05 = quantile95 = Double.NaN;
	}

	private void processValues (double[] valueSet) {
		// Remove Nans, trimming the set if needed
		double[] set = valueSet;
		java.util.Arrays.sort(set);
		int validCount = 0;
		for (int i = 0; i < set.length; i++) {
			if (Double.isNaN(set[i])) {
				break;
			}
			validCount++;
		}
		if (validCount == 0) {
			invalidate();
			return;
		}
		
		if (validCount < set.length) {
			set = Arrays.copyOf(set, validCount);
		}
		
		// Analyze Stats
		min = set[0];
		max = set[set.length - 1];

		count = set.length;
		sum = 0.0;
		for (int i = 0; i < count; i++) {
			sum += set[i];
		}
		mean = sum / count;
		
		if (count > 1) {
			double SumDSquare = 0;
			for (int i = 0; i < count; i++) {
				double d = set[i] - mean;
				SumDSquare += (d * d);
			}
			stdev = Math.sqrt(SumDSquare / (count - 1));			
			coefVar = stdev / mean;			
		} else {
			stdev = 0.0;			
			coefVar = 0.0;			
		}
		
		int medIdx1 = (set.length - 1) / 2;
		int medIdx2 = set.length / 2;
		double y1 = set[medIdx1];
		double y2 = set[medIdx2];
		median = (y1 + y2) / 2.0;
		
		int q25Idx = set.length / 4; // Skip these many at each end
		int q75Idx = (set.length - 1) - q25Idx; // Skip that many at each end
		quartile25 = set[q25Idx];
		quartile75 = set[q75Idx];
		
		int q10Idx = set.length / 10; // Skip these many at each end
		int q90Idx = (set.length - 1) - q10Idx;
		quantile10 = set[q10Idx];
		quantile90 = set[q90Idx];
		
		int q05Idx = set.length / 20; // Skip these many at each end
		int q95Idx = (set.length - 1) - q05Idx;
		quantile05 = set[q05Idx];
		quantile95 = set[q95Idx];
	}
	
	public int getCount() {
		return count;
	}

	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	public double getMean() {
		return mean;
	}

	public double getStdev() {
		return stdev;
	}
	
	public double getCoeffOfVariation() {
		return coefVar;
	}
	
	public double getSum() {
		return sum;
	}

	public double getMedian() {
		return median;
	}

	public double getQuartile25() {
		return quartile25;
	}

	public double getQuartile75() {
		return quartile75;
	}

	public double getQuantile10() {
		return quantile10;
	}

	public double getQuantile90() {
		return quantile90;
	}

	public double getQuantile05() {
		return quantile05;
	}

	public double getQuantile95() {
		return quantile95;
	}
}
