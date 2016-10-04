package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;


public class SampleReadsRetriever {
	
	private Locus[]     loci;
	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleReadsRetriever (Locus[] loci) throws AnalysisException  {
		// Parse configuration file
		this.loci = loci;
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
		getUnmappedLocusReads (sample, loci, mappedReadLists);
		
		return mappedReadLists;
	}
	
	private void getMappedLocusReads (Sample sample, Locus locus, ArrayList<MappedRead> mappedReadsList) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		GenomeRegion readSearchInterval =  locus.getReadSearchInterval();
		String chrName = readSearchInterval.getChromosome();
		String mapName = sample.getBamChromosomeMap();
		String transChrName = ChromosomeMap.getMappedChromosomeName(chrName, mapName);
		//System.out.println("Locus: " + locus.name+ " - Chr: " + chrName + " - Map: "+mapName+ " - Translated: "+transChrName);
		SAMRecordIterator it = samReader.query(transChrName, readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true);
		while (it.hasNext()) {
			SAMRecord record = it.next();
			
			// If the read is ungapped, we can use the BAM start/end coords
			if ((record.getCigarLength() == 1) && (record.getCigarString().endsWith("M"))) {
				MappedRead sr = new MappedRead(record, locus, record.getAlignmentStart(), MappedRead.MAPPED);
				mappedReadsList.add(sr);
			} else {
				// If not, take the ungapped read, and treat it as if unmapped, try to find an anchor
				@SuppressWarnings("unused")
				boolean matched = matchUnmappedReadAtLocus (record, locus, mappedReadsList, MappedRead.REMAPPED);
			}
		}
		it.close();
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
				MappedRead sr = new MappedRead(record, locus, anchors[aIdx], anchorPos, mappingStatus);
				mappedReadList.add(sr);
				return true;
		    }									
		}					
		return false;
	}
}
