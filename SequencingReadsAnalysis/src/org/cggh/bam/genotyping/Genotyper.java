package org.cggh.bam.genotyping;


public class Genotyper {
	
	protected int minCallReads;
	protected int minAlleleReads;
	protected double minAlleleProp;
	
   public Genotyper(int minCallReads, int minAlleleReads, double minAlleleProp) {
		super();
		this.minCallReads = minCallReads;
		this.minAlleleReads = minAlleleReads;
		this.minAlleleProp = minAlleleProp;
	}

	public boolean hasSufficientReads (int totalReads) {
        return (totalReads >= minCallReads);
    }
    
	public boolean isValidAllele (int alleleReads, int totalReads) {
		if (alleleReads < minAlleleReads) {
			return false;
		}
		double alleleProp = ((double) alleleReads) / ((double) totalReads);
		return (alleleProp >= minAlleleProp);
	}
}
