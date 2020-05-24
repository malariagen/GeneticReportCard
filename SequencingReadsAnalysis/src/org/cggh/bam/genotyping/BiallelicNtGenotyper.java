package org.cggh.bam.genotyping;

import org.cggh.common.counters.NtAlleleCounter;
import org.cggh.common.counters.AlleleCounter.AlleleCount;

public class BiallelicNtGenotyper {
	
	private AlleleValidator validator;
	private int minCallReads;
	
	public BiallelicNtGenotyper () {
		this (new AlleleValidator.AlleleValidatorByReadCountProportion(), 
		      GenotypingConfig.MIN_READS_FOR_CALL);		
	}
	
	public BiallelicNtGenotyper (AlleleValidator validator, int minCallReads) {
		this.validator = validator;
		this.minCallReads = minCallReads;		
	}
	
    public boolean hasSufficientReads (NtAlleleCounter ntCounts) {
    	int totalReads = ntCounts.getCumulativeCount();
        return (totalReads >= minCallReads);
    }
	
    public char getMajorityAllele (NtAlleleCounter ntCounts) {
    	int totalReads = ntCounts.getCumulativeCount();
        if (totalReads < minCallReads) {
        	return '-';
        }
	    AlleleCount[] sortedCounts = ntCounts.getSortedAlleleCounts();
	    AlleleCount top = sortedCounts[0];
        if (!validator.isValidAllele(top.getCount(), totalReads)) {
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
	public double getFrequencyCall (NtAlleleCounter ntCounts, char ref, char nref) {
	    int refReads = ntCounts.getCount(ref);
	    int nrefReads = ntCounts.getCount(nref);
	    int totalReads = refReads + nrefReads;
        if (totalReads < minCallReads) {
        	return (Double.NaN);
        }
	    if (validator.isValidAllele(refReads, totalReads)) {
		    if (validator.isValidAllele(nrefReads, totalReads)) {
		    	return ((double)nrefReads / (double)totalReads);
		    } else {
		    	return (0.0);
		    }
	    } else {
		    if (validator.isValidAllele(nrefReads, totalReads)) {
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
