package org.cggh.bam.target;

import org.cggh.common.sequence.*;
import java.util.*;

public class AminoSampleCall extends SampleCall {
	
	protected int    call;
	protected String aminoAllele;
	protected String aminoNrefAllele;
	protected String aminoAlleleSummary;
	
	// This constructor used when reading data from file.
	public AminoSampleCall (int call, String aa, String aaNref, String aaSummary, SampleCall ntCall) {
	    super (ntCall);
		this.call = call;
		this.aminoAllele = aa;
		this.aminoNrefAllele = aaNref;
		this.aminoAlleleSummary = aaSummary;
	}

	// This constructor used when computing the amino acid call from nucleotide sequence.
	public AminoSampleCall (SampleCall ntCall) {
	    super (ntCall);
	    this.call = ntCall.call;
	    if (this.call == CALL_MISSING) {
	    	return;
	    }

	    // Get the amino acids corresponding to the called codons, collapsing codons that generate the same amino
		String ntRef = ntCall.getRefSequence();
		String aaRef = SequenceUtilities.translateNtSequence(ntRef);
	    ArrayList<String> aaList = new ArrayList<String>();
	    String[] ntAlleles = allele.split(",");
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
				this.call = isRef ? CALL_WILDTYPE : CALL_MUTANT;
			} else {
				this.call = CALL_HET;
			}
			sb.append(aa);
			sbNref.append(isRef ? "." : aa);
		}
		aminoAllele = sb.toString();
		aminoNrefAllele = sbNref.toString();
		
		// Translate the sequences in the summary to amino sequences
		sb = new StringBuffer();
	    String[] summaryAlleles = alleleSummary.split(";");
		for (int aIdx = 0; aIdx < summaryAlleles.length; aIdx++) {
			if (aIdx > 0) {
				sb.append(',');
			}
		    String[] summaryParts = summaryAlleles[aIdx].split(":");
		    String ntAllele = summaryParts[0];
		    String aminoAllele = SequenceUtilities.translateNtSequence(ntAllele);
		    sb.append(aminoAllele);
		    sb.append('(');
		    sb.append(ntAllele);
		    sb.append(')');
		    sb.append(':');
		    sb.append(summaryParts[1]);
		}
		aminoAlleleSummary = sb.toString();
	}
	
	public String getAminoAllele() {
		return aminoAllele;
	}

	public String getAminoNrefAllele() {
		return aminoNrefAllele;
	}

	public String getAminoAlleleSummary() {
		return aminoAlleleSummary;
	}
	
	public static AminoSampleCall makeMissingCall () {
		return new AminoSampleCall (SampleCall.CALL_MISSING, null, null, null, SampleCall.makeMissingCall());
	}
	
}

