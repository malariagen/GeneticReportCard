package org.cggh.bam.target;

import org.cggh.common.sequence.SequenceUtilities;

public class AminoSampleCall extends SampleCall {
	
	protected String aminoAllele;
	protected String aminoNrefAllele;
	protected String aminoAlleleSummary;
	protected boolean isLenient;
	

	public AminoSampleCall (SampleCall ntCall, Target target) {
	    super (ntCall);
	    this.isLenient = ntCall.isLenient();
	    if (ntCall.call == CALL_MISSING) {
	    	return;
	    }

		String refSequence = target.getTargetRefSeq();
	    String[] ntAlleles = allele.split(",");
		StringBuffer sb = new StringBuffer();
		StringBuffer sbNref = new StringBuffer();
		for (int idx = 0; idx < ntAlleles.length; idx++) {
			if (idx > 0) {
				sb.append(',');
				sbNref.append('/');
			}
			String ntAllele = ntAlleles[idx];
			boolean isRef = refSequence.equals(ntAllele);
			String aminoAllele = SequenceUtilities.translateNtSequence(ntAllele);
			sb.append(aminoAllele);
			sbNref.append(isRef ? "." : aminoAllele);
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
		    sb.append(':');
		    sb.append(summaryParts[1]);
		}
		aminoAlleleSummary = sb.toString();
	}

	public boolean isLenient() {
		return isLenient;
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
}

