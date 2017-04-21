package org.cggh.bam;

public abstract class Genotyper {
	
	public abstract boolean isValidAllele (int alleleReads, int totalReads);
	
	public static class Genotyper5percent extends Genotyper {
		public boolean isValidAllele (int alleleReads, int totalReads) {
			// Minimum reads for the allele is 2 reads, and 5% of total
			
			// Anything below 60 reads total will require 2 reads
			if (totalReads < 60) {
				return (alleleReads >= 2);
			}
			
			// Above 60 reads total, require 5% of total
			int minReads = totalReads / 20;
			return (alleleReads >= minReads);
		}
	}
}

