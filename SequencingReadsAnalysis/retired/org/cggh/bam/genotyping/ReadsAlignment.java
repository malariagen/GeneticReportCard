package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class ReadsAlignment {
	
	public static final int MAX_DIFFERENCES_FROM_CONSENSUS = 10;

	public static final String BAM_COMMENT_PREFIX = "@CO\t";
	public static final String ALIGN_START_PREFIX = BAM_COMMENT_PREFIX+"genotype.alignment.start=";
	public static final String TIMESTAMP_PREFIX   = BAM_COMMENT_PREFIX+"genotype.timestamp=";
	public static final String TIMESTAMP_FORMAT   = "yyyy-MM-dd HH:mm";

	private DateFormat timestampFormatter = new SimpleDateFormat(TIMESTAMP_FORMAT);

	private Locus        locus;
	private Sample       sample;
	private MappedRead[] mappedReads;
	private int          alignStart;
	private int          readCount;
	private int          alignLen;
	private char[][]     alignment;
	private char[]       consensus;
	private int[]        differences;
	private int          misalignedCount;
	private Sequence     refSequence;
	
	
	public ReadsAlignment (MappedRead[] sampleReads, Locus locusConfig, Sample sample) throws AnalysisException {
		this.locus = locusConfig;
		this.mappedReads = sampleReads;
		this.sample = sample;
		this.readCount = sampleReads.length;
		//log.info("Start alignment building");

		// Align all reads according to their anchor position
		alignment = makeAlignment (sampleReads);
		//log.info("Alignment built");
		
        // Get a reference sequence as comparator
		refSequence = getReferenceSequence();
		
		// Get the consensus sequence for this alignment
		consensus = computeConsensus();
		//log.info("Computed consensus: "+new String(consensus));
		
		// Determine how good is the alignment of each read
		differences = countDifferencesFromConsensus ();
		misalignedCount = 0;
		for (int rIdx = 0; rIdx < readCount; rIdx++) {
			if (hasTooManyDifferences (rIdx)) {
				misalignedCount++;
			}
		}
	}
	
	public MappedRead[] getMappedReads() {
		return mappedReads;
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

	public int getMisalignedCount() {
		return misalignedCount;
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
	
	
	private Sequence getReferenceSequence () throws AnalysisException {
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
	
	public void outputAlignment (File outFolder) throws AnalysisException {
		Sequence[] alignSeq;
		if (mappedReads.length > 0) {
			alignSeq = new Sequence[mappedReads.length+1];
			alignSeq[0] = this.refSequence;
			for (int rIdx = 0; rIdx < mappedReads.length; rIdx++) {
				MappedRead r = mappedReads[rIdx];
				String flags = "";
				if (hasTooManyDifferences(rIdx)) {
					flags += "[MISALIGNED]";
				}
				if (!r.coversAtLeastOneTargets()) {
					flags += "[NO_TARGET_COVERAGE]";
				} else {
					if (r.hasLowQualityTargetSeq()) {
						flags += "[LOW_TARGET_Q]";
					}
				}
				if (r.isReversed()) {
					flags += "[-]";
				}
				switch(r.getMappingStatus()) {
				case MappedRead.REMAPPED:
					flags += "[R]";
					break;
				case MappedRead.UNMAPPED:
					flags += "[U]";
					break;
				}
				String seqTitle = flags.isEmpty() ? r.getId() : flags + " " + r.getId();
				alignSeq[rIdx+1] = new Sequence (seqTitle, new String(alignment[rIdx]));
			}
		} else {
			alignSeq = new Sequence[0];
		}
		
		// Write out to file the reads alignment
		File outLocusFolder = FileUtilities.checkFolder(outFolder, locus.getName(), true);
		File fastaFile = new File(outLocusFolder, sample.getName()+'-'+locus.getName()+".fasta");
		String fasta = SequenceUtilities.makeFastaAlignment(alignSeq);
		try {
			FileUtilities.writeFileContent(fasta, fastaFile);
		} catch (IOException e) {
			throw new AnalysisException ("Error writing file "+ fastaFile.getAbsolutePath()+ ": "+ e);
		}
		
		// Write out to file the SAM lines
		File samFile = new File(outLocusFolder, sample.getName()+'-'+locus.getName()+".sam");
		BufferedTextOutput samOut = new BufferedTextOutput(new UncompressedOutputTextStore (samFile), 16 * 1024 * 1024);
		StringBuffer sb = samOut.getStringBuffer();
		sb.append(TIMESTAMP_PREFIX+timestampFormatter.format(new Date()));
		sb.append('\n');
		sb.append(ALIGN_START_PREFIX+alignStart);
		sb.append('\n');
		sb.append(BAM_COMMENT_PREFIX+"Flag\tRefName\tPos\tMapQ\tCigar\tRefNext\tPosNext\tTmplLen\tSeq\tQual");
		for (int rIdx = 0; rIdx < mappedReads.length; rIdx++) {
			MappedRead sr = mappedReads[rIdx];
			sb.append('\n');
			sb.append(sr.getSamString());
			samOut.commitIfBufferFull();
		}
		samOut.close();
	}
}
