package org.cggh.bam.genotyping;

import org.cggh.bam.*;

public class Genotyper {
	
	protected int minCallReads;
	protected int minAlleleReads;
	protected double minAlleleProp;
	
	public Genotyper(BaseAnalysisConfig config) {
	    this.minCallReads   = config.getMinCallReadCount();
		this.minAlleleReads = config.getMinAlleleReadCount();
		this.minAlleleProp  = config.getMinAlleleReadProp();
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
