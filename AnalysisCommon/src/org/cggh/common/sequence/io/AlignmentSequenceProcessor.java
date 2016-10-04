package org.cggh.common.sequence.io;

public class AlignmentSequenceProcessor implements SequenceProcessor {

	private int lineLength = -1;
	
	/**
	 * Prepares an aligned sequence after it is read.
	 * It ensures all sequences are of equal length, and uses whitespace to pad the 
	 * leading and trailing positions
	 * 
	 * @param seq The sequence being prepared
	 * @return a copy of the string after preparation
	 */
	public String processSequence(String seq) throws SequenceSourceException {
		
		// Check the length
		if (lineLength < 0) {
			lineLength = seq.length();
		} else if (lineLength != seq.length()) {
			throw new SequenceSourceException (
					"Inconsistent sequence (length is "+seq.length() +" instead of "+lineLength+")");
		}
		
		
		// Find the start of sequence data
		int pos;
		for (pos = 0; pos < seq.length(); pos++) {
			if (isSequenceChar (seq.charAt(pos))) {
				break;
			}
		}
		int seqStart = pos;
		
		// Find the end of sequence data
		for (pos = (seq.length() -1); pos > seqStart; pos--) {
			if (isSequenceChar (seq.charAt(pos))) {
				break;
			}
		}
		int seqEnd = pos;

		// Now re-write the sequence, using spaces for leading and trailing
		StringBuffer sb = new StringBuffer(seq.length());		
		for (int i = 0; i < seq.length(); i++) {
			if ((i < seqStart) || (i > seqEnd)) {
				sb.append(' ');
			} else {
				sb.append(seq.charAt(i));
			}
		}
		return sb.toString();
	}
	
	private boolean isSequenceChar (char c) {
		switch(c) {
		case ' ':
		case '-':
		case '*':
		case 'X':
		case 'x':
			return false;
		}
		return true;
	}
	
}
