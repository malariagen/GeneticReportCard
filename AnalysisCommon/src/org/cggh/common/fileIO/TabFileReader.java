package org.cggh.common.fileIO;

import java.io.*;

public class TabFileReader extends DelimitedReader {
	public TabFileReader (Reader reader) throws IOException {
		this(reader, DEFAULT_BUFFER_SIZE);
	}
	
	public TabFileReader (Reader reader, int bufferSize) throws IOException {
		super(reader, bufferSize);
		setDelimiter ("\t");
	}
}