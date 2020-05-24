package org.cggh.bam.genotyping;

public abstract class AlleleValidator {
	
	public abstract boolean isValidAllele (int alleleReads, int totalReads);
	
	public static class AlleleValidatorByReadCountProportion extends AlleleValidator {
		
		protected double readsPropMultiplier;
		protected int    totalReadsThreshold;
		
	    public AlleleValidatorByReadCountProportion () {
			this(GenotypingConfig.ALLELE_READS_PROPORTION_THRESHOLD);
		}
		
		public AlleleValidatorByReadCountProportion (double proportionThreshold) {
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

}

