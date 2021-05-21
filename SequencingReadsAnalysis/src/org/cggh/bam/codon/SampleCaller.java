package org.cggh.bam.codon;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.sequence.SequenceUtilities;

import java.util.*;


public class SampleCaller {
	
	private int     minCallReadCount;
	private int     minAlleleReadCount;
	private double  minAlleleReadProp;
	
	public SampleCaller(BaseAnalysisConfig config) {
		this.minCallReadCount   = config.getMinCallReadCount();
		this.minAlleleReadCount = config.getMinAlleleReadCount();
		this.minAlleleReadProp  = config.getMinAlleleReadProp();
	}

	public SampleCall callSample (AlignmentTarget target, LabelCounters ntAlleleCounters) {
		String ref = target.getTargetRefSeq();
		
		ArrayList<String> validAlleleList = new ArrayList<String>();
		int totalReads = ntAlleleCounters.getTotal();
		if (totalReads >= minCallReadCount) {
			// We have enough reads to make a call (potentially).
			// However, we need to discard alleles that do not have enough coverage and therefore we suspect could be noise.
			// Start with the allele read counts, ordered by read count (descending)
			int totalValidReads = 0;
		    LabelCounter[] counters = ntAlleleCounters.getSortedCounters();
			for (int cIdx = 0; cIdx < counters.length; cIdx++) {
				// In order for an allele to be accepted, we need a minimum of reads, and these must make up at least a minimum proportion of the total reads.
				boolean validAllele = isValidAllele (counters[cIdx], totalReads);
				if (validAllele) {
					validAlleleList.add(counters[cIdx].getLabel());
					totalValidReads += counters[cIdx].getCount();
				} else {
					// We've found an allele that does not meet criteria; since the array is sorted by read counts, any subsequent allele will also not be valid.
					// So we terminate allele checking. Before we exit, we need to reassess whether we have enough good reads to make a call, after eliminating
					// all the counters that we think might have been noise. If we do not have enough, reset the called alleles before exiting.
					if (totalValidReads < minCallReadCount) {
						validAlleleList.clear();
					}
					break;
				}
			}
		}
		
		String alleleSummary = getAlleleSummary (ntAlleleCounters);
		SampleCall ntCall = new SampleCall (SampleCall.CALL_MISSING, ref, null, null, alleleSummary);
		int alleleCount = validAlleleList.size();
		if (alleleCount == 1) {
			String allele = validAlleleList.get(0);
			boolean isRef = ref.equals(allele);
			int call = isRef ? SampleCall.CALL_WILDTYPE : SampleCall.CALL_MUTANT;
			String nrefAllele = isRef ? "." : allele;
			ntCall = new SampleCall (call, ref, allele, nrefAllele, alleleSummary);
		} else if (alleleCount > 1) {
			StringBuffer sb = new StringBuffer();
			StringBuffer sbNref = new StringBuffer();
			for (int cIdx = 0; cIdx < alleleCount; cIdx++) {
				if (cIdx > 0) {
					sb.append(',');
					sbNref.append('/');
				}
				String a = validAlleleList.get(cIdx);
				boolean isRef = ref.equals(a);
				sb.append(a);
				sbNref.append(isRef ? "." : a);
			}
			String allele = sb.toString();
			String nrefAllele = sbNref.toString();
			ntCall = new SampleCall (SampleCall.CALL_HET, ref, allele, nrefAllele, alleleSummary);
		}
		return ntCall;
	}
	
	// Translate the nt sequences in the summary to amino sequences, and incorporate them into the summary
	private String getAlleleSummary (LabelCounters ntAlleleCounters) {
		LabelCounter[] counters = ntAlleleCounters.getSortedCounters();
		if (counters.length == 0) {
			return "-";
		}
		
		StringBuffer sb = new StringBuffer();
		for (int aIdx = 0; aIdx < counters.length; aIdx++) {
			if (aIdx > 0) {
				sb.append(',');
			}
		    String ntAllele = counters[aIdx].getLabel();
		    int count = counters[aIdx].getCount();
		    String aminoAllele = SequenceUtilities.translateNtSequence(ntAllele);
		    sb.append(aminoAllele);
		    sb.append('(');
		    sb.append(ntAllele);
		    sb.append(')');
		    sb.append(':');
		    sb.append(count);
		}
		return sb.toString();
	}
	
	//
	// In order for an allele to be accepted, we need a minimum of reads, and these must make up at least a minimum proportion of the total reads.
    //
	private boolean isValidAllele (LabelCounter counter, int totalReads) {
		int alleleReads = counter.getCount();
		if (alleleReads < minAlleleReadCount) {
			return false;
		}
		double alleleProp = ((double) alleleReads) / ((double) totalReads);
		return (alleleProp >= minAlleleReadProp);
	}
}
