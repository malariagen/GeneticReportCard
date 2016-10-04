package org.cggh.bam.genotyping;

import org.cggh.common.counters.*;


public class SampleCall {
	
	public static final int CALL_MISSING  = 1;
	public static final int CALL_WILDTYPE = 2;
	public static final int CALL_MUTANT   = 3;
	public static final int CALL_HET      = 4;
	
	int call;
	String aaAllele;
	String aaNrefAllele;
	String aaAlleleSummary;
	
	public SampleCall (LabelCounters aCounters, Target t, int minCallReads, int minAlleleCallReads) {
		
		aaAlleleSummary = aCounters.getSummary();
		LabelCounter[] counters = aCounters.getSortedCounters();
		
		boolean missing = (counters.length == 0) || (counters[0].getCount() < minAlleleCallReads) || (aCounters.getTotal() < minCallReads);
		if (missing) {
			call = CALL_MISSING;
			return;
		}
		
		if (counters.length == 1) {
			String aa = ((TargetAlleleCounter)counters[0]).aaSequence;
			boolean isRef = t.getTargetAaSeq().equals(aa);
			call = isRef ? CALL_WILDTYPE : CALL_MUTANT;
			aaAllele = aa;
			aaNrefAllele = isRef ? "." : aa;
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
			String aa = ((TargetAlleleCounter)counters[cIdx]).aaSequence;
			boolean isRef = t.getTargetAaSeq().equals(aa);
			sb.append(aa);
			sbNref.append(isRef ? "." : aa);
		}
		aaAllele = sb.toString();
		aaNrefAllele = sbNref.toString();
	}
	
	public boolean isMissing() {
		return (call == CALL_MISSING);
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
}

