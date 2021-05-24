package org.cggh.bam.genotyping;

import org.cggh.common.counters.*;

public class LabelCounterGenotyper extends Genotyper {
	
    public LabelCounterGenotyper(int minCallReads, int minAlleleReads, double minAlleleProp) {
		super(minCallReads, minAlleleReads, minAlleleProp);
	}
    
	public boolean hasSufficientReads (LabelCounters counts) {
        return hasSufficientReads(counts.getTotal());
    }
    
    public String getMajorityAllele (LabelCounters counts) {
        if (!hasSufficientReads(counts)) {
        	return "-";
        }
    	int totalReads = counts.getTotal();
    	LabelCounter[] counters = counts.getSortedCounters();
    	LabelCounter top = counters[0];
        if (!isValidAllele(top.getCount(), totalReads)) {
        	return "-";
        }
	    return top.getLabel();
   }
}
