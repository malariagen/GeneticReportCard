package org.cggh.common.sequence.io;

public interface SequenceProcessor {

	/**
	 * Prepares a sequernce after it is read.
	 * 
	 * @param seq The sequence being prepared
	 * @return a copy of the string after preparation (may be the string itself if no preparation needed)
	 */
	public String processSequence(String seq) throws SequenceSourceException;
	
}
