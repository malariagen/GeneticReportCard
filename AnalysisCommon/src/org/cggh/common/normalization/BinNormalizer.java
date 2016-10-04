package org.cggh.common.normalization;

import java.io.*;

public class BinNormalizer {

	private NormalizerBin[]  bins;
	
//	public BinningNormalizer (double minX, double maxX, int binCount) {
//		double[] binMaxValues = new double[binCount];
//		double binWidth = (maxX - minX)/(double)binCount;
//		double binMax = minX;
//		for (int i = 0; i < binCount; i++) {
//			binMax += binWidth;
//			binMaxValues[i] = binMax;
//		}
//		bins = initializeBins (binMaxValues) ;
//	}
	
	public BinNormalizer (double[] binMaxValues) {
		bins = initializeBins (binMaxValues) ;
	}
	
	private NormalizerBin[] initializeBins (double[] binMaxValues) {
		NormalizerBin[] newBins = new NormalizerBin[binMaxValues.length];
		for (int i = 0; i < binMaxValues.length; i++) {
			newBins[i] = new NormalizerBin(binMaxValues[i]);
		}
		return newBins;
	}
	
	public void addValue(BinNormalizable value) {
		if (value != null) {
			NormalizerBin bin = findBin(value);
			bin.addValue(value);
		}
	}
	
	public NormalizerBin findBin(BinNormalizable value) {
		double binningValue = value.getBinningValue();
		for (int i = 0; i < bins.length; i++) {
			if (binningValue < bins[i].getMaxX()) {
				return bins[i];
			}
		}
		return bins[bins.length-1];
	}
	
	public void normalize() {
		for (int i = 0; i < bins.length; i++) {
			bins[i].normalize();
		}
	}
	
	public NormalizerBin[] getBins() {
		return bins;
	}

	public void serialize (Writer w) {
		PrintWriter pw = new PrintWriter(w);
		pw.print("Num\tBinMax\tMean\tStdev\tCount");
		for (int i = 0; i < bins.length; i++) {
			NormalizerBin bin = bins[i];
			pw.print("\n" + (i+1)
				   + "\t" +bin.getMaxX()
				   + "\t" +bin.getMeanY()
				   + "\t" +bin.getStdevY()
				   + "\t" +bin.getValueCount()
				   );
		}
		pw.close();
	}
	

	public static class NormalizerBin extends Normalizer {
		private double maxX;
		
		public NormalizerBin (double maxX) {
			super();
			this.maxX = maxX;
		}

		public double getMaxX() {
			return maxX;
		}
	}
}
