package org.cggh.common.subsampling;

public class BinCounters {
	private double[] maxBinValues;
	private int[]    counts;
	
	public BinCounters (double[] maxBinValues) {
		this.maxBinValues = maxBinValues;
		this.counts = new int[maxBinValues.length];
	}
	
	public int getNumBins () {
		return counts.length;
	}
	
	public int getBinCount (int binIdx) {
		return counts[binIdx];
	}
	
	public int[] getBinCounts () {
		return counts;
	}
	
	public double getBinMaxValue (int binIdx) {
		return maxBinValues[binIdx];
	}
	
	public double[] getBinMaxValues () {
		return maxBinValues;
	}
	
	public void incrementBinForValue (double value) {
		int binIdx = getBinIdx (value);
		if (binIdx >= 0) {
			counts[binIdx]++;
		}
	}
	
	public void incrementBinByIndex (int binIdx) {
		if (binIdx >= 0) {
			counts[binIdx]++;
		}
	}
	
	public void reset () {
		for (int idx = 0; idx < counts.length; idx++) {
			counts[idx] = 0;
		}
	}
	
	public int getBinIdx (double value) {
		if (value > 0.0) {
			for (int idx = 0; idx < maxBinValues.length; idx++) {
				if (value <= maxBinValues[idx]) {
					return idx;
				}
			}
			//if (value == maxBinValues[maxBinValues.length-1]) {
			//	return (maxBinValues.length-1);
			//}
		}
		return -1;
	}
}
