package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.sequence.*;
import java.util.*;

public class ReadsAlignment {
	
	public static final int MAX_DIFFERENCES_FROM_CONSENSUS = 10;

	private Locus        locus;
	private Sample       sample;
	private MappedRead[] alignedReads;
	private MappedRead[] misalignedReads;
	private int          alignStart;
	private int          readCount;
	private int          alignLen;
	private char[][]     alignment;
	private char[]       consensus;
	private int[]        differences;
	private Sequence     referenceSequence;
	
	
	public ReadsAlignment (Sample sample, Locus locus, MappedRead[] sampleReads) throws AnalysisException {
		this.sample = sample;
		this.locus = locus;
		//log.info("Start alignment building");

		// Align all reads according to their anchor position
		readCount = sampleReads.length;
		alignment = makeAlignment (sampleReads);
		//log.info("Alignment built");
		
        // Get a reference sequence as comparator
		referenceSequence = retrieveReferenceSequence();
		
		// Get the consensus sequence for this alignment
		consensus = computeConsensus();
		//log.info("Computed consensus: "+new String(consensus));
		
		// Determine how good is the alignment of each read
		differences = countDifferencesFromConsensus ();
		
		ArrayList<MappedRead> alignedReadList = new ArrayList<MappedRead>(sampleReads.length);
		ArrayList<MappedRead> misalignedReadList = new ArrayList<MappedRead>();
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			if (hasTooManyDifferences (rIdx)) {
				misalignedReadList.add(sampleReads[rIdx]);
			} else {
				alignedReadList.add(sampleReads[rIdx]);
			}
		}
		alignedReads = alignedReadList.toArray(new MappedRead[alignedReadList.size()]);
		misalignedReads = misalignedReadList.toArray(new MappedRead[misalignedReadList.size()]);
	}
	
	public Sample getSample() {
		return sample;
	}

	public MappedRead[] getAllMappedReads() {
		MappedRead[] sampleReads = new MappedRead[alignedReads.length+misalignedReads.length];
		System.arraycopy(alignedReads, 0, sampleReads, 0, alignedReads.length);
		System.arraycopy(misalignedReads, 0, sampleReads, alignedReads.length, misalignedReads.length);
		return sampleReads;
	}

	public MappedRead[] getAlignedReads() {
		return alignedReads;
	}

	public MappedRead[] getMisalignedReads() {
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
		return (differences[rIdx] > MAX_DIFFERENCES_FROM_CONSENSUS);
	}

	private char[][] makeAlignment (MappedRead[] sampleReads) {
		// Compute the maximum extent of the alignemnt
		int alignEnd = alignStart = -1;
		for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
			MappedRead sr = sampleReads[rIdx];
			int rStart = sr.getStartPos();
			int rEnd = rStart + sr.getSequence().length() - 1;
			if (rIdx == 0) {
				alignStart = rStart;
				alignEnd = rEnd;
			} else {
				alignStart = (alignStart < rStart) ? alignStart : rStart;
				alignEnd = (alignEnd > rEnd) ? alignEnd : rEnd;
			}
		}
		alignLen = 1 + alignEnd - alignStart;

		// Create an array for storing the alignments
		char[][] result = new char[readCount][alignLen];
		
		// Write the read sequences in the right places
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			Arrays.fill(result[rIdx], '-');
			MappedRead sr = sampleReads[rIdx];
			int idx = sr.getStartPos() - alignStart;
			String seq = sr.getSequence();
			for (int i = 0; i < seq.length(); i++) {
				result[rIdx][idx+i] = seq.charAt(i);
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
