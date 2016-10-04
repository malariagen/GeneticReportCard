package org.cggh.common.sequence.io;

import java.io.*;

import org.cggh.common.sequence.*;

public abstract class SequenceReader implements FileTypeProcessor {

	public static final Class<?>[] FILE_PROCESSOR_CLASSES = {
		TabbedSequenceReader.class,
		//RawSequenceReader.class,
		FastaSequenceReader.class,
	};
	
	/**
	 * Reads the next sequence from a LineNumberReader. This will be implemented by
	 * subclasses in a format-specific way.
	 * 
	 * @return a Sequence object, retrieved and prepared for analysis, or null if none could be found
	 * @throws SequenceSourceException if there are file operation errors
	 */
	public abstract Sequence readNextSequence (LineNumberReader reader) throws SequenceSourceException;
		
	
	/**
	 * Gets the next line of text from the file, skipping lines that 
	 * are blank or only whitespace
	 * 
	 * @return the next valid line, or null if EOF
	 * @throws IOException if there is an IO Exception in the low-level file reading
	 */
	protected String getNextValidLine(LineNumberReader reader) throws SequenceSourceException {
		String s = null;
		while (true) {
			
			// Advance to next line		
			try {
				s = reader.readLine();
			} catch (IOException e) {
				throw new SequenceSourceException("Error reading valid line", e, reader.getLineNumber());
			}
			
			// End of file: return null
			if (s == null) {
				break;
			}
			
			// Check the line has at least some non-whitespace characters, else skip
			if ((s.length() != 0) && (s.trim().length() != 0)) {
				break;
			}
		}
		return s;
	}	
		
}
