package org.cggh.bam;

import org.cggh.bam.target.TargetAnalysisConfig;
import org.cggh.common.exceptions.*;
import org.cggh.common.sequence.*;
import java.util.*;

public class ReadsAlignment {
	
	private static int maxReadMismatches = TargetAnalysisConfig.DEFAULT_MAX_READ_MISMATCHES;
	
	public static void configure (LocusAnalysisConfig config) {
		ReadsAlignment.maxReadMismatches  = config.getMaxReadMismatches();	
	}
	
	private Locus    locus;
	private Sample   sample;
	private Read[]   alignedReads;
	private Read[]   misalignedReads;
	private int      alignStart;
	private int      readCount;
	private int      alignLen;
	private char[][] alignment;
	private char[]   consensus;
	private int[]    differences;
	private Sequence referenceSequence;

	public ReadsAlignment (Sample sample, Locus locus, Read[] sampleReads) throws AnalysisException {
		this.sample = sample;
		this.locus = locus;
		//log.info("Start alignment building");

		// Align all reads according to their anchor position
		readCount = sampleReads.length;
		if (readCount == 0) {
			return;
		}
		alignment = makeAlignment (sampleReads, locus);
		//log.info("Alignment built");
		
        // Get a reference sequence as comparator
		referenceSequence = retrieveReferenceSequence();
		
		// Get the consensus sequence for this alignment
		consensus = computeConsensus();
		//log.info("Computed consensus: "+new String(consensus));
		
		// Determine how good is the alignment of each read
		differences = countDifferencesFromConsensus ();
		
		ArrayList<Read> alignedReadList = new ArrayList<Read>(sampleReads.length);
		ArrayList<Read> misalignedReadList = new ArrayList<Read>();
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			if (hasTooManyDifferences (rIdx)) {
				misalignedReadList.add(sampleReads[rIdx]);
			} else {
				alignedReadList.add(sampleReads[rIdx]);
			}
		}
		alignedReads = alignedReadList.toArray(new Read[alignedReadList.size()]);
		misalignedReads = misalignedReadList.toArray(new Read[misalignedReadList.size()]);
	}
	
	public Sample getSample() {
		return sample;
	}

	public Read[] getAllMappedReads() {
		Read[] sampleReads = new Read[alignedReads.length+misalignedReads.length];
		System.arraycopy(alignedReads, 0, sampleReads, 0, alignedReads.length);
		System.arraycopy(misalignedReads, 0, sampleReads, alignedReads.length, misalignedReads.length);
		return sampleReads;
	}

	public Read[] getAlignedReads() {
		return alignedReads;
	}

	public Read[] getMisalignedReads() {
		return misalignedReads;
	}

	public Sequence getReferenceSequence() {
		return referenceSequence;
	}

	public int getAlignStart() {
		return alignStart;
	}

	public int getReadCount() {
		return readCount;
	}

	public int getAlignLen() {
		return alignLen;
	}

	public String getAlignedReadSequence(int rIdx) {
		return new String(alignment[rIdx]);
	}

	private int[] countDifferencesFromConsensus () {
		int[] diffs = new int[readCount];
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			diffs[rIdx] = 0;
			for (int i = 0; i < alignLen; i++) {
				char nt = alignment[rIdx][i];
				if ((nt != '-') && (nt != consensus[i])) {
					diffs[rIdx]++;
				}
			}
		}
		return diffs;
	}
	
	public boolean hasTooManyDifferences (int rIdx) {
		return (differences[rIdx] > maxReadMismatches);
	}

	private char[][] makeAlignment (Read[] sampleReads, Locus locus) {
		// Compute the maximum extent of the alignemnt
		int alignEnd = alignStart = Integer.MIN_VALUE;
		for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
			Read sr = sampleReads[rIdx];
			int rStart = sr.getStartPos();
			int rEnd = rStart + sr.getSequence().length() - 1;
			if (rIdx == 0) {
				alignStart = rStart;
				alignEnd = rEnd;
			} else {
				alignStart = (alignStart < rStart) ? alignStart : rStart;
				//if (rStart < 1) { System.out.println("Here"); } // This may happen sometimes in amplicons
				alignEnd = (alignEnd > rEnd) ? alignEnd : rEnd;
			}
		}
		
		// Trim the alignment to fit in the locus region being investigated
		int locusStart = locus.getReadSearchInterval().getStartPos();
		if (alignStart < locusStart) {
			alignStart = locusStart;
		}
		int locusEnd = locus.getReadSearchInterval().getStopPos();
		if (alignEnd > locusEnd) {
			alignEnd = locusEnd;
		}
		
		// Create an array for storing the alignments+
		alignLen = 1 + alignEnd - alignStart;
		char[][] result = new char[readCount][alignLen];
		
		// Write the read sequences in the right places
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			Arrays.fill(result[rIdx], '-');
			Read sr = sampleReads[rIdx];
			String seq = sr.getSequence();
			for (int i = 0; i < seq.length(); i++) {
				int pos = (sr.getStartPos() + i);
				if ((pos < alignStart) || (pos > alignEnd)) {
					continue;
				}
				int posIdx = (pos - alignStart);
				result[rIdx][posIdx] = seq.charAt(i);
			}
		}
		return result;
	}
	
	private Sequence retrieveReferenceSequence () throws AnalysisException {
		if (readCount == 0) {
			return null;
		}
		String chrName = locus.getReadSearchInterval().getChromosome();
		Sequence chrSeq = ReferenceGenome.getChrSequence(chrName);
		String refSeq = chrSeq.getData().substring(alignStart-1, alignStart+alignLen-1);
		String refSeqTitle = "REF|"+chrName+":"+alignStart+"-"+(alignStart+alignLen-1);
		return new Sequence(refSeqTitle, refSeq);
	}
	
	
	private char[] computeConsensus () {
		char[] result = new char[alignLen];
		int[] ntCounts = new int[4];
		for (int i = 0; i < alignLen; i++) {
			ntCounts[0] = ntCounts[1] = ntCounts[2] = ntCounts[3] = 0;
			for (int j = 0; j < readCount; j++) {
				char nt = alignment[j][i];
				if ((nt != '-') && (nt != 'N')) {
					int ntIdx = nt2Index (nt);
					ntCounts[ntIdx]++;
				}
			}
			result[i] = getMajorityAllele (ntCounts);
		}
		return result;
	}
	
	private char getMajorityAllele (int[] ntCounts) {
		return (index2Nt(getMajorityIndex(ntCounts)));
	}
	
	private int getMajorityIndex (int[] ntCounts) {
		int mCount = -1;
		int mIdx = -1;
		for (int i = 0; i < 4; i++) {
			if (ntCounts[i] > mCount) {
				mIdx = i;
				mCount = ntCounts[i];
			}
		}
		return mIdx;
	}
	
	private int nt2Index (char nt) {
		switch (nt) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		}
		return -1;
	}
	private char index2Nt (int index) {
		switch (index) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		}
		return (char)0;
	}
}
