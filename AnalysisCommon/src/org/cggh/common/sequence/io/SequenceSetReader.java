package org.cggh.common.sequence.io;

import java.io.*;
import java.util.*;

import org.cggh.common.sequence.Sequence;


/**
 * Inputs a set of sequences. No post processing of the sequences is performed.
 * Will read in files containing sequences in different formats, depending on the SequenceReader selected.
 * 
 * @author olivo
 */

public class SequenceSetReader {
	
	private SequenceReader seqReader;
	private SequenceProcessor processor;
	
	
	public SequenceSetReader (SequenceReader seqReader) {
		this(seqReader, new NullSequenceProcessor());
	}
	
	protected SequenceSetReader (SequenceReader seqReader, SequenceProcessor processor) {
		this.processor = processor;
		this.seqReader = seqReader;
	}
	
	

	/**
	 * Reads a stream and return the aligned sequences as strings, after they are
	 * prepared for entropy analysis.
	 * 
	 * @param inFile the input file, in the format specified
	 * @return an array containing the sequences retrieved and prepared for analysis
	 * @throws IOException
	 */
	public Sequence[] readSequences (InputStream inStream) throws SequenceSourceException {

		// Process the file according to format
		ArrayList<Sequence> v = new ArrayList<Sequence>();
		LineNumberReader reader = new LineNumberReader(new InputStreamReader(inStream));
		reader.setLineNumber(1);
		try {
			int count = 0;
			while (true) {

				// Get sequence data line
				Sequence sequence = seqReader.readNextSequence (reader);
				if (sequence == null) {
					break;
				}

				String seqData = sequence.getData();
				String seqId = sequence.getId();
				
				seqId = (seqId != null) ? seqId : "S"+(count++);
				seqData = processor.processSequence(seqData);
				
				// Add the sequence to the list
				v.add(new Sequence(seqId, seqData));
			}
		} catch (Exception e) {
			throw new SequenceSourceException("Error reading sequences", e, reader.getLineNumber());
			
		} finally {
			try {
				reader.close();
			} catch (IOException e2) {
				throw new SequenceSourceException("Error closing input file", e2);
			}
		}
 
		// Package the results into an array of strings
		Sequence[] sequences = new Sequence[v.size()];
		sequences = (Sequence[]) v.toArray(sequences);

		return sequences;
	}

	/**
	 * Reads a file and return the aligned sequences as strings, after they are
	 * prepared for entropy analysis.
	 * 
	 * @param inFile the input file, in the format specified
	 * @return an array containing the sequences retrieved and prepared for analysis
	 * @throws IOException
	 */
	public Sequence[] readSequences (File inFile) throws SequenceSourceException {
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(inFile);
		} catch (IOException e) {
			throw new SequenceSourceException("Error opening input file", e);
		}
		return readSequences(fis);
	}


	
}