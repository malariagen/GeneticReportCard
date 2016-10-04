package org.cggh.bam.genotyping;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.cggh.common.exceptions.*;
import org.cggh.common.sequence.*;
import org.cggh.common.sequence.io.*;
import org.cggh.common.util.*;

import htsjdk.samtools.*;
import java.io.*;
import java.text.*;
import java.util.*;


public class ReadsAlignmentRetriever {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());

	private SimpleDateFormat timestampFormatter = new SimpleDateFormat(ReadsAlignment.TIMESTAMP_FORMAT);
	
	private Locus   locusConfig;
	private File          locusFolder;
	
	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

	public ReadsAlignmentRetriever (Locus locusConfig, File lociFolder) throws AnalysisException {
		this.locusConfig = locusConfig;
		this.locusFolder = FileUtilities.checkFolder(lociFolder, locusConfig.getName(), false);
	}

	public File getFastaFile (Sample sample) throws AnalysisException {
		String filenameRoot =  sample.getName()+"-"+locusConfig.getName();
		File fastaFile = new File (locusFolder, filenameRoot+".fasta");
		return fastaFile;
	}
	
	public File getSamFile (Sample sample) throws AnalysisException {
		String filenameRoot =  sample.getName()+"-"+locusConfig.getName();
		File samFile = new File (locusFolder, filenameRoot+".sam");
		return samFile;
	}
	
	public ReadsAlignment retrieveReadsAlignment (Sample sample) throws AnalysisException {

		File samFile = getSamFile (sample);
		File fastaFile = getFastaFile (sample);

		int alignStart = -1;
		@SuppressWarnings("unused")
		Date timestamp = null;
		
		Sequence[] readAlignSeqs = null;
		try {
			SequenceSetReader seqSetReader = new SequenceSetReader (new FastaSequenceReader());
			readAlignSeqs = seqSetReader.readSequences(fastaFile);
		} catch (Exception e) {
			log.error ("Error reading reads alignment FASTA file for sample "+sample.getName()+": "+e);
			return null;
		}
		if (readAlignSeqs.length == 0) {
			// Empty alignment
			return null;
		}
		// Remove REF sequence
		if (readAlignSeqs[0].getId().startsWith("REF|")) {
			Sequence[] cleanReadAlignSeqs = new Sequence[readAlignSeqs.length-1];
			System.arraycopy(readAlignSeqs, 1, cleanReadAlignSeqs, 0, cleanReadAlignSeqs.length);
			readAlignSeqs = cleanReadAlignSeqs;
		}
		
		MappedRead[] mappedReads = new MappedRead[readAlignSeqs.length];
		try {
			SamReader samReader = samReaderFactory.open(samFile);
			
			SAMFileHeader header = samReader.getFileHeader();
			List<String> comments = header.getComments();
			for (String comment : comments) {
				if (comment.startsWith(ReadsAlignment.ALIGN_START_PREFIX)) {
					String alignStartStr = comment.substring(ReadsAlignment.ALIGN_START_PREFIX.length(), comment.length());
					alignStart = Integer.parseInt(alignStartStr.trim());
				}
				if (comment.startsWith(ReadsAlignment.TIMESTAMP_PREFIX)) {
					String timestampStr = comment.substring(ReadsAlignment.TIMESTAMP_PREFIX.length(), comment.length());
					timestamp = timestampFormatter.parse(timestampStr.trim());
				}
			}
			
			int rIdx = 0;
			SAMRecordIterator it = samReader.iterator();
			while (it.hasNext()) {
				SAMRecord record = it.next();
				
				String seqId = readAlignSeqs[rIdx].getId();
				int mappingStatus = MappedRead.MAPPED;
				if (seqId.contains("[R]")) {
					mappingStatus = MappedRead.REMAPPED;
				} else if (seqId.contains("[U]")) {
					mappingStatus = MappedRead.UNMAPPED;
				}
				
				String aSeq = readAlignSeqs[rIdx].getData();
				int readStartPos = alignStart;
				for (int i = 0; i < aSeq.length(); i++) {
					if (aSeq.charAt(i) != '-') {
						readStartPos += i;
						break;
					}
				}
				mappedReads[rIdx] = new MappedRead(record, locusConfig, readStartPos, mappingStatus);
				rIdx++;
			}
			it.close();
			samReader.close();
			
		} catch (Exception e) {
			throw new AnalysisException ("Error parsing SAM file headers at locus "+ locusConfig.getName()+ " for sample "+sample.getName()+": "+e);
		}
		if (alignStart < 0) {
			throw new AnalysisException ("No alignment start information in SAM file headers at locus "+ locusConfig.getName()+ " for sample "+sample.getName());
		}
		
		ReadsAlignment ra = new ReadsAlignment (mappedReads, locusConfig, sample);
		return ra;
	}
}
