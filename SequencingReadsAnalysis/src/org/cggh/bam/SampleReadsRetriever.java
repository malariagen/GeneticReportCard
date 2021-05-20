package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;


public class SampleReadsRetriever {
	
	private Locus[]           loci;
	private boolean           skipUnmappedReadsAnalysis;
	private int               maxIndelSize;

	private SamReaderFactory  samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleReadsRetriever (BamConfig config) throws AnalysisException  {
		// Parse configuration file
		this.loci = config.getLoci();
		this.skipUnmappedReadsAnalysis = config.getSkipUnmappedReadsAnalysis();
		this.maxIndelSize = config.getMaxIndelSize();
	}

	@SuppressWarnings("unchecked")
	public ArrayList<MappedRead>[] retrieveSampleReads (Sample sample) throws AnalysisException, IOException  {
		// Read the reads from the SAM file
		ArrayList<MappedRead>[] mappedReadLists = new ArrayList[loci.length];
		
		for (int i = 0; i < loci.length; i++) {
			Locus locus = loci[i];
			mappedReadLists[i] = new ArrayList<MappedRead>();
			getMappedLocusReads (sample, locus, mappedReadLists[i]);
		}
		
		// Then search unmapped reads (this does not necessarily find them all, but will do)
		if (!skipUnmappedReadsAnalysis) {
			getUnmappedLocusReads (sample, loci, mappedReadLists);
		}
		
		return mappedReadLists;
	}
	
	private void getMappedLocusReads (Sample sample, Locus locus, ArrayList<MappedRead> mappedReadsList) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		GenomeRegion readSearchInterval =  locus.getReadSearchInterval();
		String chrName = readSearchInterval.getChromosome();

		SAMRecordIterator it = samReader.query(chrName, readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true);
		while (it.hasNext()) {
			SAMRecord record = it.next();
			// This will ignore all secondary and supplementary alignments, reads not passing vendor filters, as well as any duplicates
			if (record.getFlags() >= 256) {
				continue;
			}
			
			// Do an initial mapping of the read
			MappedRead sr = new MappedRead(record, locus, record.getAlignmentStart());

			// Unless the read is mapped "as is" (e.g. CIGAR string is something like "150M"), process the CIGAR to refine mapping against the reference
			if ((record.getCigarLength() > 1) || (!record.getCigarString().endsWith("M"))) {
				try {
					applyCigar (sr, record.getCigar());
				} catch (CigarException e) {
					sr.unmap();
				}
			}
			if (sr.getMappingStatus() == MappedRead.MAPPED) {
				mappedReadsList.add(sr);
			} else if (!skipUnmappedReadsAnalysis) {
				// If not, take the ungapped read, and treat it as if unmapped, try to find an anchor
				@SuppressWarnings("unused")
				boolean matched = matchUnmappedReadAtLocus (record, locus, mappedReadsList, MappedRead.REMAPPED);
			}
		}
		it.close();
	}

	private static class CigarException extends AnalysisException {
	    public CigarException (String msg) {
	        super (msg);
	    }
	}
	
	private void applyCigar(MappedRead sr, Cigar cigar) throws CigarException {
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

	private void getUnmappedLocusReads (Sample sample, Locus[] loci, ArrayList<MappedRead>[] mappedReadLists) throws AnalysisException {
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

	private boolean matchUnmappedRead (SAMRecord record, Locus[] loci, ArrayList<MappedRead>[] mappedReadLists) throws AnalysisException {
		boolean matched = false;
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus locus = loci[lIdx];
			if (!locus.getAnalyzeUnmappedReads()) {
				continue;
			}
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			if (matchUnmappedReadAtLocus (record, locus, mappedReadList, MappedRead.UNMAPPED)) {
				matched = true;  // The same unmapped read may have anchors that match multiple loci, so do not give up after finding a match
			}
		}
		return matched;
	}

	private boolean matchUnmappedReadAtLocus (SAMRecord record, Locus locus, ArrayList<MappedRead> mappedReadList, int mappingStatus) throws AnalysisException {
		// Does the read contain an anchor?
		String readSequence = record.getReadString();
		Anchor[] anchors = locus.getAnchors();
		for (int aIdx = 0; aIdx < anchors.length; aIdx++) {
			Matcher m = anchors[aIdx].getRegex().matcher(readSequence);
		    if (m.find()) {
		    	int anchorPos = m.start();
				MappedRead sr = new MappedRead(record, locus, anchors[aIdx], anchorPos);
				sr.setMappingStatus (mappingStatus);
				mappedReadList.add(sr);
				return true;
		    }									
		}					
		return false;
	}
}
