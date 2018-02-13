package org.cggh.common.fileIO;

import java.io.*;

public class DelimitedReader extends LineNumberReader {
	
    public static int DEFAULT_BUFFER_SIZE = 8192;
    
    public static String DEFAULT_DELIMITER = ",";
    public static String DEFAULT_COMMENT_PREFIX = "#";
	
	private String   delimiter     = DEFAULT_DELIMITER;
	private String   commentPrefix = DEFAULT_COMMENT_PREFIX;
	private String   currLine = null;
	private String[] currFields = null;
	
	public DelimitedReader (Reader reader) throws IOException {
		this(reader, DEFAULT_BUFFER_SIZE);
	}
	
	public DelimitedReader (Reader reader, int bufferSize) throws IOException {
		super(reader, bufferSize);
		this.setLineNumber(1);
		this.currLine = null;
	}
	
	public String getCurrentLine () {
		return currLine;
	}
	
	public String[] getCurrentLineFields () {
		return currFields;
	}
	
	public void setDelimiter (String delimiter) {
		this.delimiter = delimiter;
	}
	
	public void setCommentPrefix (String commentPrefix) {
		this.commentPrefix = commentPrefix;
	}
	
	public String[] getNextValidLine() throws IOException {
		currFields = null;
		boolean discardComments = (commentPrefix != null);
		while (true) {
			
			// Advance to next line
			currLine = readLine();
			
			// End of file: return null
			if (currLine == null) {
				break;
			}
			
			// Check the line has at least some non-whitespace characters, else skip
			if ((currLine.length() != 0) && (currLine.trim().length() != 0)) {
				if (discardComments && currLine.startsWith(commentPrefix)) {
					continue;
				}
				currFields = currLine.split(delimiter);
				break;
			}
		}
		return currFields;
	}
}