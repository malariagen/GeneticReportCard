package org.cggh.common.textStore;

import java.io.*;


public abstract class OutputTextStore extends TextStore {
	
	public OutputTextStore (File folder, String filename) {
		super (folder, filename);
	}
	
	public Writer getWriter() throws IOException {
		return getWriter(false);
	}
		
	public abstract Writer getWriter(boolean append) throws IOException;
	
	public abstract void closeWriter() throws IOException ;
}
