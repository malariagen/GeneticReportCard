package org.cggh.common.sequence.io;


import java.io.IOException;
import java.io.LineNumberReader;

import org.cggh.common.sequence.Sequence;

/**
 *  Processor for Tab delimited sequences, i.e. one line = name<tab>sequence
 */
public class TabbedSequenceReader extends SequenceReader {
	
	public String getDescription () {
		return "ID<tab>Sequence";
	}
	
	public String[] getFileExtensions () {
		return new String[]{"taln"};
	}
	
	/**
	 * Reads the next sequence from a file of sequences in a tab delimited format: 
	 * one line = name<tab>sequence
	 * Skips blank lines.
	 * 
	 * @return a sequence string, retrieved and prepared for analysis
	 * @throws IOException
	 */
	public Sequence readNextSequence (LineNumberReader reader) throws SequenceSourceException {

		// Get and check FASTA header line
		String line = getNextValidLine(reader);
		if (line == null) {
			return null;
		} 
		
		// Get and check sequence data line (second part of the entry)
		String[] parts = line.split("\t");
		if (parts.length != 2) {
			throw new SequenceSourceException("Line does not contain exactly one tab", reader.getLineNumber());
		}
		return new Sequence (parts[0].trim(), parts[1].trim());
	}
}

