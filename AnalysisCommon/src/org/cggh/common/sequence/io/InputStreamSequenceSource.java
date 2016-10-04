package org.cggh.common.sequence.io;


import java.io.*;

import org.cggh.common.sequence.Sequence;

public class InputStreamSequenceSource extends SequenceSource {
	
	private Sequence[] sequences;
	
	public InputStreamSequenceSource (InputStream inStream, SequenceSetReader reader) throws SequenceSourceException {
		try {
			sequences = reader.readSequences(inStream);							
		} finally {
			try {
				inStream.close();
			} catch (IOException e) {}
		}
	}
	
	public Sequence[] getSequences() throws SequenceSourceException {
		return sequences;
	}
	
	public String getName() {
		return "<Unknown Input Stream>";
	}
}
