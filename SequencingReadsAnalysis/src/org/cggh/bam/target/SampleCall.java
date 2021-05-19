package org.cggh.bam.target;


public class SampleCall {
	
	public static final int CALL_MISSING  = 1;
	public static final int CALL_WILDTYPE = 2;
	public static final int CALL_MUTANT   = 3;
	public static final int CALL_HET      = 4;
	
	protected int    call;
	protected String ref;
	protected String allele;
	protected String nrefAllele;
	protected String alleleSummary;
	
	/*
	 * ONLY USED BY SUPERCLASSES (AminoSampleCall)
	 */
	protected SampleCall (SampleCall ntCall) {
		this (ntCall.call, ntCall.ref, ntCall.allele, ntCall.nrefAllele, ntCall.alleleSummary);
	}
	
	public SampleCall (int call, String ref, String allele, String nrefAllele, String alleleSummary) {
		this.call = call;
		this.allele = allele;
		this.nrefAllele = nrefAllele;
		this.alleleSummary = alleleSummary;
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

	public String getRefSequence() {
		return allele;
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
	
	public static int getCallFromString(String callStr) {
		if (callStr.equals("WT")) {
			return CALL_WILDTYPE;
		} else if (callStr.equals("MI")) {
			return CALL_MISSING;
		} else if (callStr.equals("MU")) {
			return CALL_MUTANT;
		} else if (callStr.equals("HE")) {
			return CALL_HET;
		}
		return -1;
	}

	public static SampleCall makeMissingCall () {
		return new SampleCall (CALL_MISSING, null, null, null, null);
	}
}