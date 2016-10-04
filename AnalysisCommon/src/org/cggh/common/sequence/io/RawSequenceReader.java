package org.cggh.common.sequence.io;


import java.io.IOException;
import java.io.LineNumberReader;

import org.cggh.common.sequence.Sequence;

/**
 *  Processor for Raw sequences, i.e. one line = one aligned sequence
 */
public class RawSequenceReader extends SequenceReader {
	
	public String getDescription () {
		return "Raw Sequences";
	}
	
	public String[] getFileExtensions () {
		return new String[]{"txt"};
	}
	
	/**
	 * Simply returns the next non-empty valid line, assuming it is a sequence
	 * Skips blank lines.
	 * 
	 * @return a sequence string, retrieved and prepared for analysis
	 * @throws IOException
	 */
	public Sequence readNextSequence (LineNumberReader reader) throws SequenceSourceException {
		String line = getNextValidLine(reader);
		if (line == null) {
			return null;
		}
		return new Sequence(null, line);  // No id
	}
}


