package org.cggh.bam.genotyping;

import org.cggh.bam.*;
import org.cggh.common.counters.AlleleCounter;
import org.cggh.common.counters.AlleleCounter.AlleleCount;

public class NucleotideGenotyper extends Genotyper {
	
    public NucleotideGenotyper(BaseAnalysisConfig config) {
		super(config);
	}

	public boolean hasSufficientReads (AlleleCounter aCounts) {
        return hasSufficientReads(aCounts.getCumulativeCount());
    }
    
    public char getMajorityAllele (AlleleCounter aCounts) {
        if (!hasSufficientReads(aCounts)) {
        	return '-';
        }
    	int totalReads = aCounts.getCumulativeCount();
	    AlleleCount[] sortedCounts = aCounts.getSortedAlleleCounts();
	    AlleleCount top = sortedCounts[0];
        if (!isValidAllele(top.getCount(), totalReads)) {
        	return '-';
        }
	    return top.getAllele();
    }
    
   	/**
   	 * Note: this ignores 3rd/4th alleles! It assumes that they are handled elsewhere
   	 * @param ntCounts
   	 * @param ref
   	 * @param nref
   	 * @return
   	 */
   	public double getFrequencyCall (AlleleCounter ntCounts, char ref, char nref) {
   	    int refReads = ntCounts.getCount(ref);
   	    int nrefReads = ntCounts.getCount(nref);
   	    int totalReads = refReads + nrefReads;
           if (totalReads < minCallReads) {
           	return (Double.NaN);
           }
   	    if (isValidAllele(refReads, totalReads)) {
   		    if (isValidAllele(nrefReads, totalReads)) {
   		    	return ((double)nrefReads / (double)totalReads);
   		    } else {
   		    	return (0.0);
   		    }
   	    } else {
   		    if (isValidAllele(nrefReads, totalReads)) {
   		    	return (1.0);
   		    }	    	
   	    }
   	    return (Double.NaN); // We should never be here...
   	}
   	
   	public int frequencyToMajorityAlleleNum (double freq) {
   		int num = frequencyToAlleleNum (freq);
   		if (num == 3) {
   			num = (freq > 0.5) ? 2 : 1;
   		}
   		return num;
   	}

   	public int frequencyToAlleleNum (double freq) {
   	    if (Double.isNaN(freq)) {
   	    	return 0;
   	    }
   	    if (freq == 0.0) return 1;
   	    if (freq == 1.0) return 2;
   	    return 3;  // Het		
   	}
}
