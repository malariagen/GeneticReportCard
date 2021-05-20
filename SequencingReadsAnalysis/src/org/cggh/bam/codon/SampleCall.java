package org.cggh.bam.codon;

import java.util.ArrayList;

import org.cggh.common.sequence.SequenceUtilities;

public class SampleCall {
	
	public static final int CALL_MISSING  = 1;
	public static final int CALL_WILDTYPE = 2;
	public static final int CALL_MUTANT   = 3;
	public static final int CALL_HET      = 4;
	
	private int    call;
	private String ref;
	
	private String ntAllele;
	private String ntNrefAllele;
	private String aaAllele;
	private String aaNrefAllele;
	private String alleleSummary;
	
	public SampleCall (int call, String ref, String ntAllele, String ntNrefAllele, String alleleSummary) {
		this.call = call;
		this.ref = ref;
		
		this.ntAllele = ntAllele;
		this.ntNrefAllele = ntNrefAllele;
		this.alleleSummary = alleleSummary;
		
		if (call != CALL_MISSING) {
			callAminoAcid ();
		}
		
		// Translate the sequences in the summary to amino sequences
	}
	
	protected SampleCall (SampleCall srcCall) {
		this.call = srcCall.call;
		this.ref = srcCall.ref;
		
		this.ntAllele = srcCall.ntAllele;
		this.ntNrefAllele = srcCall.ntNrefAllele;
		
		this.aaAllele = srcCall.aaAllele;
		this.aaNrefAllele = srcCall.aaNrefAllele;
		
		this.alleleSummary = srcCall.alleleSummary;
	}
	
	private void callAminoAcid () {
	    // Get the amino acids corresponding to the called codons, collapsing codons that generate the same amino
		String ntRef = getRefSequence();
		String aaRef = SequenceUtilities.translateNtSequence(ntRef);
	    ArrayList<String> aaList = new ArrayList<String>();
	    String[] ntAlleles = ntAllele.split(",");
		for (int idx = 0; idx < ntAlleles.length; idx++) {
			String nt = ntAlleles[idx];
			String aa = SequenceUtilities.translateNtSequence(nt);
			if (!aaList.contains(aa)) {
				aaList.add(aa);
			}
		}
	    
	    // Create strings showing the calles amino(s), also in "nonref" notation
		StringBuffer sb = new StringBuffer();
		StringBuffer sbNref = new StringBuffer();
		for (int idx = 0; idx < aaList.size(); idx++) {
			String aa = aaList.get(idx);
			boolean isRef = aaRef.equals(aa);
			if (idx > 0) {
				sb.append(',');
				sbNref.append('/');
				this.call = CALL_HET;
			} else {
				this.call = isRef ? CALL_WILDTYPE : CALL_MUTANT;
			}
			sb.append(aa);
			sbNref.append(isRef ? "." : aa);
		}
		aaAllele = sb.toString();
		aaNrefAllele = sbNref.toString();
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
		return ntAllele;
	}

	public String getNtAllele() {
		return ntAllele;
	}

	public String getNtNrefAllele() {
		return ntNrefAllele;
	}

	public String getAminoAllele() {
		return aaAllele;
	}

	public String getAminoNrefAllele() {
		return aaNrefAllele;
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
		return new SampleCall (SampleCall.CALL_MISSING, null, null, null, null);
	}
	
}