package org.cggh.common.sequence.io;

/**
 * Inputs a set of sequences. No post processing of the sequences is performed.
 * Will read in files containing sequences in different formats, depending on the SequenceReader selected.
 * 
 * @author olivo
 */

public class AlignmentReader extends SequenceSetReader {
	
	public AlignmentReader (SequenceReader seqReader) {
		super(seqReader, new AlignmentSequenceProcessor());
	}
	
}