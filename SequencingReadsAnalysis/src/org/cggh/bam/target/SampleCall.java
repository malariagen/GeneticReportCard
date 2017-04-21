package org.cggh.bam.target;

import org.cggh.bam.Genotyper;
import org.cggh.common.counters.*;

public class SampleCall {
	
	public static final int CALL_MISSING  = 1;
	public static final int CALL_WILDTYPE = 2;
	public static final int CALL_MUTANT   = 3;
	public static final int CALL_HET      = 4;
	
	protected int    call;
	protected String allele;
	protected String nrefAllele;
	protected String alleleSummary;
	protected Genotyper genotyper = new Genotyper.Genotyper5percent();
	
	public SampleCall (LabelCounters aCounters, Target t, int minCallReads) {
		
		alleleSummary = aCounters.getSummary();
		LabelCounter[] counters = aCounters.getSortedCounters();

		boolean missing = true;
		int totalReads = aCounters.getTotal();
		int alleleCount = 0;
		if (totalReads >= minCallReads) {
			for (int cIdx = 0; cIdx < counters.length; cIdx++) {
				boolean valid = genotyper.isValidAllele(counters[cIdx].getCount(), totalReads);
				if (!valid) {
					break;
				}
				missing = false;
				alleleCount++;
			}
		}
		if (missing) {
			call = CALL_MISSING;
			return;
		}
		
		if (alleleCount == 1) {
			allele = counters[0].getLabel();
			boolean isRef = t.getTargetRefSeq().equals(allele);
			call = isRef ? CALL_WILDTYPE : CALL_MUTANT;
			nrefAllele = isRef ? "." : allele;
			return;
		}
		
		call = CALL_HET;
		StringBuffer sb = new StringBuffer();
		StringBuffer sbNref = new StringBuffer();
		for (int cIdx = 0; cIdx < alleleCount; cIdx++) {
			if (cIdx > 0) {
				sb.append(',');
				sbNref.append('/');
			}
			String a = counters[cIdx].getLabel();
			boolean isRef = t.getTargetRefSeq().equals(a);
			sb.append(a);
			sbNref.append(isRef ? "." : a);
		}
		allele = sb.toString();
		nrefAllele = sbNref.toString();
	}
	
	/*
	 * ONLY USED BY SUPERCLASSES (AminoSampleCall)
	 */
	protected SampleCall (SampleCall ntCall) {
		this.call = ntCall.call;
		this.allele = ntCall.allele;
		this.nrefAllele = ntCall.nrefAllele;
		this.alleleSummary = ntCall.alleleSummary;
	}
	
	public boolean isMissing() {
		return (call == CALL_MISSING);
	}
	
	public boolean isLenient() {
		return false;
	}
	
	public int getCall() {
		return call;
	}

	public String getCallString() {
		switch (call) {
		case CALL_MISSING:  return "MI";
		case CALL_WILDTYPE: return "WT";
		case CALL_MUTANT:   return "MU";
		case CALL_HET:      return "HE";
		}
		return null;
	}

	public String getAllele() {
		return allele;
	}

	public String getNrefAllele() {
		return nrefAllele;
	}

	public String getAlleleSummary() {
		return alleleSummary;
	}
	
	
	public static class LenientSampleCall extends SampleCall {
		public LenientSampleCall (LabelCounters aCounters, Target t, int minCallReads) {
			super (aCounters, t, minCallReads);
		}
		public boolean isLenient() {
			return true;
		}
	}
}