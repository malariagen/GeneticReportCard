package org.cggh.common.sequence.io;

import java.io.IOException;
import java.io.LineNumberReader;

import org.cggh.common.sequence.Sequence;

/**
 *  Processor for FASTA sequences, i.e. one fasta header, followed by 
 *  one or more line containing the raw aligned sequence
 */
public class FastaSequenceReader extends SequenceReader {
	
	private String lastHeader = null;

	public String getDescription () {
		return "FASTA sequences";
	}
	
	public String[] getFileExtensions () {
		return new String[]{"afa", "fa","fas","fasta"};
	}
	
	/**
	 * Reads the next sequence from a file of sequences in FASTA format: 
	 * one fasta header, followed by one or more lines containing the 
	 * raw aligned sequence. Skips blank lines.
	 * 
	 * @return a sequence string, retrieved and prepared for analysis
	 * @throws IOException
	 */
	public Sequence readNextSequence (LineNumberReader reader) throws SequenceSourceException {
		
		// If it's the first record (or we're after the last one)
		// then there is no previous header. So we need to find 
		// the first header, or find out there is no next record.
		// 
		if (lastHeader == null) {
			lastHeader = findFastaHeader (reader);
			if (lastHeader == null) {
				return null;
			}
		}
		String seqId = lastHeader.substring(1); // Remove '>'
		
		// We've got a header- go through the sequence data 
		// until we find another header
		StringBuffer sb = new StringBuffer();
		String nextLine = getNextValidLine(reader);
		while ((nextLine != null) && !isFastaHeader(nextLine)) {
			sb.append(nextLine.trim());
			nextLine = getNextValidLine(reader);
		}
		lastHeader = nextLine;
		String seqData = sb.toString();
		if (seqData.length() == 0) {
			throw new SequenceSourceException("Missing sequence data", reader.getLineNumber());
		}
		return new Sequence(seqId, seqData);
	}
	
	private String findFastaHeader (LineNumberReader reader) throws SequenceSourceException {
		// Get and check FASTA header line
		String nextLine = getNextValidLine(reader);
		if (nextLine == null) {
			return null;
		}
		if (!isFastaHeader(nextLine)) {
			throw new SequenceSourceException("Invalid FASTA Header", reader.getLineNumber());
		}
		return nextLine;
	}
	
	private boolean isFastaHeader (String line) {
		return (line.trim().startsWith(">"));
	}
}
