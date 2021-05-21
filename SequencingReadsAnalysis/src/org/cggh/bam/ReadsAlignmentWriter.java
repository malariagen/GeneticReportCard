package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class ReadsAlignmentWriter {
	
	public static final String BAM_COMMENT_PREFIX = "@CO\t";
	public static final String ALIGN_START_PREFIX = BAM_COMMENT_PREFIX+"genotype.alignment.start=";
	public static final String TIMESTAMP_PREFIX   = BAM_COMMENT_PREFIX+"genotype.timestamp=";
	public static final String TIMESTAMP_FORMAT   = "yyyy-MM-dd HH:mm";

	private DateFormat timestampFormatter = new SimpleDateFormat(TIMESTAMP_FORMAT);


	public void outputAlignment (ReadsAlignment ra, Locus locus, File outFolder) throws AnalysisException {
		Sequence[] alignSeq;
		Read[] mappedReads = ra.getAllMappedReads();
		if (mappedReads.length > 0) {
			alignSeq = new Sequence[mappedReads.length+1];
			alignSeq[0] = ra.getReferenceSequence();
			for (int rIdx = 0; rIdx < mappedReads.length; rIdx++) {
				Read r = mappedReads[rIdx];
				String flags = "";
				if (ra.hasTooManyDifferences(rIdx)) {
					flags += "[MISALIGNED]";
				}
				if (r.isReversed()) {
					flags += "[-]";
				}
				switch(r.getMappingStatus()) {
				case Read.ANCHORED:
					flags += "[A]";
					break;
				case Read.UNMAPPED:
					flags += "[U]";
					break;
				}
				String seqTitle = flags.isEmpty() ? r.getId() : flags + " " + r.getId();
				alignSeq[rIdx+1] = new Sequence (seqTitle, ra.getAlignedReadSequence(rIdx));
			}
		} else {
			alignSeq = new Sequence[0];
		}
		
		// Write out to file the reads alignment
		Sample sample = ra.getSample();
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
		sb.append(ALIGN_START_PREFIX+ra.getAlignStart());
		sb.append('\n');
		sb.append(BAM_COMMENT_PREFIX+"Flag\tRefName\tPos\tMapQ\tCigar\tRefNext\tPosNext\tTmplLen\tSeq\tQual");
		for (int rIdx = 0; rIdx < mappedReads.length; rIdx++) {
			Read sr = mappedReads[rIdx];
			sb.append('\n');
			sb.append(sr.getSamString());
			samOut.commitIfBufferFull();
		}
		samOut.close();
	}
}
