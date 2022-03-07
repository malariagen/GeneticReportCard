package org.cggh.bam.genotyping;

import org.cggh.bam.BaseAnalysisConfig;

public class AlleleValidator {
	
	protected double readsPropMultiplier;
	protected int    totalReadsThreshold;
	
	public AlleleValidator (BaseAnalysisConfig config) {
		double proportionThreshold = config.getMinAlleleReadProp();
		
		// This is a rounding trick to make sure the threshold is included in the valid interval.
		// e.g. if threshold is 5%, and total reads are 60, then min reads should be 3, not 4 
		// (which would be the case without this tiny adjustment)
		readsPropMultiplier = proportionThreshold - 1E-6;
		
		// Min number of total reads; below this value, the min reads for a call is 2
		totalReadsThreshold = (int)(2.0 / proportionThreshold);
	}
	
	public boolean isValidAllele (int alleleReads, int totalReads) {
		// Anything below minTotalReads reads total will require 2 reads to call an allele
		int minAlleleReads = (totalReads <= totalReadsThreshold) ? 2 :  (1 + (int)(((double)totalReads) * readsPropMultiplier));
		return (alleleReads >= minAlleleReads);
	}
}

