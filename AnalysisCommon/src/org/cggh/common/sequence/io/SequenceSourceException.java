package org.cggh.common.sequence.io;

/**
 * Exception class for this alignmennt processor.
 */
public class SequenceSourceException extends Exception {

	public SequenceSourceException() {
	}
	
	public SequenceSourceException(String msg) {
		super(msg);
	}
	
	public SequenceSourceException(String msg, int lineNumber) {
		super(buildMessage(msg, lineNumber));
	}
	
	public SequenceSourceException(String msg, Throwable cause) {
		super(msg, cause);
	}
	
	public SequenceSourceException(String msg, Throwable cause, int lineNumber) {
		super(buildMessage(msg, lineNumber), cause);
	}
	
	private static String buildMessage(String msg, int lineNumber) {
		StringBuffer sb = new StringBuffer();
		if (lineNumber >= 0) {
			sb.append("Line ");
			sb.append(lineNumber);
			sb.append(": ");
		}
		sb.append(msg);
		return sb.toString();
	}
}
