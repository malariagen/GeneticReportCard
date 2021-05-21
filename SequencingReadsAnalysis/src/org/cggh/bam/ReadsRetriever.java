package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;


public class ReadsRetriever {
	
	private Locus[]           loci;
	private boolean           analyzeUnmappedReads;
	private int               maxIndelSize;
	private boolean           useAlignment;

	private SamReaderFactory  samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public ReadsRetriever (BaseAnalysisConfig config) throws AnalysisException  {
		// Parse configuration file
		this.loci = config.getLoci();
		this.useAlignment = config.getUseBamAlignment();
		this.analyzeUnmappedReads = config.getAnalyzeUnmappedReads();
		this.maxIndelSize = config.getMaxIndelSize();
		
		// Verify we only have single search intervals for alignmet-based tasks
		if (useAlignment) {
			for (int lIdx = 0; lIdx < loci.length; lIdx++) {
				Locus locus = loci[lIdx];
				GenomeRegion[] searchIntervals = locus.getReadSearchIntervals();
				if (searchIntervals.length > 1) {
					throw new AnalysisException ("Error at locus "+locus.getName()+": an alignment-base task cannot use more than one search interval per locus.");
				}
			}
		}
	}

	@SuppressWarnings("unchecked")
	public ArrayList<Read>[] retrieveSampleReads (Sample sample) throws AnalysisException, IOException  {
		// Read the reads from the SAM file
		ArrayList<Read>[] readLists = new ArrayList[loci.length];
		
		for (int i = 0; i < loci.length; i++) {
			Locus locus = loci[i];
			readLists[i] = new ArrayList<Read>();
			GenomeRegion[] searchIntervals = locus.getReadSearchIntervals();
			for (int j = 0; j < searchIntervals.length; j++) {
				GenomeRegion interval = searchIntervals[j];
				getMappedLocusReads (sample, locus, interval, readLists[i], useAlignment);
			}
		}
		
		// Then search unmapped reads (this does not necessarily find them all, but will do)
		if (analyzeUnmappedReads) {
			getUnmappedLocusReads (sample, loci, readLists);
		}
		
		return readLists;
	}
	
	private void getMappedLocusReads (Sample sample, Locus locus, GenomeRegion readSearchInterval, ArrayList<Read> readsList, boolean useAlignment) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		String chrName = readSearchInterval.getChromosome();

		SAMRecordIterator it = samReader.query(chrName, readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true);
		while (it.hasNext()) {
			SAMRecord record = it.next();
			// This will ignore all secondary and supplementary alignments, reads not passing vendor filters, as well as any duplicates
			if (record.getFlags() >= 256) {
				continue;
			}
			
			if (useAlignment) { 
				// Use the BAM alignment to do an initial mapping of the read
				Read sr = Read.createMappedRead(record, locus, record.getAlignmentStart());
				
				// Unless the read is mapped "as is" (e.g. CIGAR string is something like "150M"), 
				// process the CIGAR to refine mapping against the reference
				if ((record.getCigarLength() > 1) || (!record.getCigarString().endsWith("M"))) {
					try {
						applyCigar (sr, record.getCigar());
					} catch (CigarException e) {
						sr.unmap();
					}
				}
				
				// If the mapping is still valid after applying CIGAR, use the alignment.
				// If not, take the ungapped read, and treat it as if unmapped, try to find an anchor
				if (sr.getMappingStatus() == Read.MAPPED) {
					readsList.add(sr);
				} else if (analyzeUnmappedReads) {
					@SuppressWarnings("unused")
					boolean matched = matchReadAtLocus (record, locus, readsList);
				}
			} else {
				// Do not use the BAM alignment, just find an anchor.
				@SuppressWarnings("unused")
				boolean matched = matchReadAtLocus (record, locus, readsList);
			}
		}
		it.close();
	}

	private static class CigarException extends AnalysisException {
	    public CigarException (String msg) {
	        super (msg);
	    }
	}
	
	private void applyCigar(Read sr, Cigar cigar) throws CigarException {
		String sequence = sr.getSequence();
		String quality = sr.getQuality();
		
		StringBuffer sequenceSb = new StringBuffer(sequence.length());
		StringBuffer qualitySb = new StringBuffer(quality.length());
		
		int seqPos = 0;
		List<CigarElement> ceList = cigar.getCigarElements();
		for (CigarElement ce : ceList) {
			CigarOperator op = ce.getOperator();
			int len = ce.getLength();
			if (op.isIndel() && (len > maxIndelSize)) {
				throw new CigarException ("Found large indel: "+len+op.name());
			}
			if (op.isAlignment()) {							// M, =, X
				// Keep all matched positions
				sequenceSb.append(sequence.substring(seqPos, seqPos+len));
				qualitySb.append(quality.substring(seqPos, seqPos+len));
				seqPos += len;
			} else if (op == CigarOperator.INSERTION) {		// I
				// Skip insertions, since they do not map against the reference
				seqPos += len;
			} else if (op == CigarOperator.DELETION) {		// D
				// Insert gaps for deletions
				for (int i = 0; i < len; i++) {
					sequenceSb.append('-');
					qualitySb.append('0');
				}
			} else if (op == CigarOperator.SOFT_CLIP) {		// S
				// Skip soft clips, start position should not change
				seqPos += len;
			} else if (op == CigarOperator.HARD_CLIP) {		// H
				// Ignore hard clips, they are not in the sequence
			} else if (op == CigarOperator.SKIPPED_REGION) { // N
				// This really should not happen!
				throw new CigarException ("Found skipped region in CIGAR string- could not process");
			} else if (op == CigarOperator.PADDING) {		// P
				// Insert gaps for padding, though it should not happen at all
				for (int i = 0; i < len; i++) {
					sequenceSb.append('-');
					qualitySb.append(' ');
				}
			} else {
				throw new CigarException ("Unkonwn element in CIGAR string: "+len+op.name());
			}
		}
		sr.updateSequence(sequenceSb.toString(), qualitySb.toString());
	}

	private void getUnmappedLocusReads (Sample sample, Locus[] loci, ArrayList<Read>[] mappedReadLists) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		SAMRecordIterator it = samReader.queryUnmapped();
		while (it.hasNext()) {
			SAMRecord record = it.next();
			boolean matched = matchUnmappedRead (record, loci, mappedReadLists);
			if (!matched) {
				SAMRecordUtil.reverseComplement(record);
				matched = matchUnmappedRead (record, loci, mappedReadLists);
			}
		}
		it.close();
	}

	private boolean matchUnmappedRead (SAMRecord record, Locus[] loci, ArrayList<Read>[] mappedReadLists) throws AnalysisException {
		boolean matched = false;
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus locus = loci[lIdx];
			if (!locus.getAnalyzeUnmappedReads()) {
				continue;
			}
			ArrayList<Read> mappedReadList = mappedReadLists[lIdx];
			if (matchReadAtLocus (record, locus, mappedReadList)) {
				matched = true;  // The same unmapped read may have anchors that match multiple loci, so do not give up after finding a match
			}
		}
		return matched;
	}

	private boolean matchReadAtLocus (SAMRecord record, Locus locus, ArrayList<Read> mappedReadList) throws AnalysisException {
		// Does the read contain an anchor?
		String readSequence = record.getReadString();
		Anchor[] anchors = locus.getAnchors();
		for (int aIdx = 0; aIdx < anchors.length; aIdx++) {
			Matcher m = anchors[aIdx].getRegex().matcher(readSequence);
		    if (m.find()) {
		    	int anchorPos = m.start();
				Read sr = Read.createAnchoredRead (record, locus, anchors[aIdx], anchorPos);
				mappedReadList.add(sr);
				return true;
		    }									
		}					
		return false;
	}
}
