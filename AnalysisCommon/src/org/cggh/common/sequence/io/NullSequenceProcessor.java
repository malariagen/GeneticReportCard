package org.cggh.common.sequence.io;

public class NullSequenceProcessor implements SequenceProcessor {

	/**
	 * Prepares a sequence after it is read, when the string requires no further processing
	 * 
	 * @param seq The sequence being prepared
	 * @return a copy of the string after preparation (may be the string itself if no preparation needed)
	 */
	public String processSequence(String seq) {
		return seq;
	}
	
}
