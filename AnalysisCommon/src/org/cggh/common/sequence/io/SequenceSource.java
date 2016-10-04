package org.cggh.common.sequence.io;

import org.cggh.common.sequence.Sequence;

public abstract class SequenceSource {
	
	public abstract Sequence[] getSequences() throws SequenceSourceException;
	
	public abstract String getName();

}
