package org.cggh.common.sequence.io;


import java.io.*;

import org.cggh.common.sequence.Sequence;

public class FileSequenceSource extends SequenceSource {
	
	private File alignmentFile;
	private SequenceSetReader reader;
	private Sequence[] sequences;
	
	public FileSequenceSource (File alignmentFile, SequenceSetReader reader) {
		this.alignmentFile = alignmentFile;
		this.reader = reader;
	}
	
	public FileSequenceSource (String alignmentFilename, SequenceSetReader reader) {
		this(new File(alignmentFilename), reader);
	}
	
	public Sequence[] getSequences() throws SequenceSourceException {
		if (sequences == null) {
		    FileInputStream inStream = null;
			try {
			    inStream = new FileInputStream(alignmentFile);
				InputStreamSequenceSource src = new InputStreamSequenceSource(inStream, reader);
				sequences = src.getSequences();
			} catch (SequenceSourceException e) {
				throw e;
			} catch (Exception e) {
				throw new SequenceSourceException (
						"Error loading alignment file "+alignmentFile.getAbsolutePath(), e);
			} finally {
				if (inStream != null) {
					try {
						inStream.close();
					} catch (Exception e2) {
						throw new SequenceSourceException (
								"Error closing alignment file "+alignmentFile.getAbsolutePath(), e2);
					}
				}
			}
		}
		return sequences;
	}
	
	public void resetSequences() {
		sequences = null;
	}

	public File getAlignmentFile() {
		return alignmentFile;
	}

	public String getName() {
		return alignmentFile.getAbsolutePath();
	}

}
