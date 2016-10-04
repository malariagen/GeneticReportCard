package org.cggh.bam.target;

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
	
	public SampleCall (LabelCounters aCounters, Target t, int minCallReads, int minAlleleCallReads) {
		
		alleleSummary = aCounters.getSummary();
		LabelCounter[] counters = aCounters.getSortedCounters();
		
		boolean missing = (counters.length == 0) || (counters[0].getCount() < minAlleleCallReads) || (aCounters.getTotal() < minCallReads);
		if (missing) {
			call = CALL_MISSING;
			return;
		}
		
		if (counters.length == 1) {
			allele = counters[0].getLabel();
			boolean isRef = t.getTargetRefSeq().equals(allele);
			call = isRef ? CALL_WILDTYPE : CALL_MUTANT;
			nrefAllele = isRef ? "." : allele;
			return;
		}
		
		call = CALL_HET;
		StringBuffer sb = new StringBuffer();
		StringBuffer sbNref = new StringBuffer();
		for (int cIdx = 0; cIdx < counters.length; cIdx++) {
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
	 * ONLY USED BY SUPERCLASSES (AminoSamnpleCall)
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
		public LenientSampleCall (LabelCounters aCounters, Target t, int minCallReads, int minAlleleCallReads) {
			super (aCounters, t, minCallReads, minAlleleCallReads);
		}
		public boolean isLenient() {
			return true;
		}
	}
}