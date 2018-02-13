package org.cggh.common.fileIO;

import java.io.*;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.InputTextStore;

public class TableInput {
	
	private InputTextStore inTextStore;
	private TabFileReader r = null;
	private String[] fieldNames;
	
	public TableInput (InputTextStore inTextStore) throws AnalysisException {
		this (inTextStore, "\t", DelimitedReader.DEFAULT_COMMENT_PREFIX);
	}
	
	public TableInput (InputTextStore inTextStore, String delimiter, String commentPrefix) throws AnalysisException {
		
		// Create a Tab-separated text file line-by-line reader
		this.inTextStore = inTextStore;
		try {
			r = new TabFileReader(inTextStore.getReader());
			r.setDelimiter(delimiter);
			r.setCommentPrefix(commentPrefix);
		} catch (IOException e) {
			throw new AnalysisException("Error opening datafile " + inTextStore.getPath() +	": "+e);
		}
		
		// Read the fields from the next line
		fieldNames = getNextValidLine();
		if (fieldNames == null) {
			throw new AnalysisException("Error opening datafile " + inTextStore.getPath() +	": no headers found");
		}
	}
	
	public TableInput (File inFile) throws AnalysisException {
		this(new InputTextStore (inFile.getParentFile(), inFile.getName()));
	}
	
	public TableInput (File inFile, String delimiter, String commentPrefix) throws AnalysisException {
		this(new InputTextStore (inFile.getParentFile(), inFile.getName()), delimiter, commentPrefix);
	}
	
	public String[] getFieldNames () {
		return fieldNames;
	}
	
	public int getFieldIndex (String fieldName) {
		for (int i = 0; i < fieldNames.length; i++) {
			if (fieldName.equals(fieldNames[i])) {
				return i;
			}
		}
		return -1;
	}
	
	public String[] getNextValidLine() throws AnalysisException {
		try {
			return r.getNextValidLine();
		} catch (IOException e) {
			throw new AnalysisException("Error reading record from file " + inTextStore.getPath() + " at line " + r.getLineNumber() + ": "+e);				
		}			
	}
	
	public int getLineNumber () {
		return r.getLineNumber();
	}
	
	public String getCurrentLine () {
		return r.getCurrentLine();
	}
	
	public String[] getCurrentLineFields () {
		return r.getCurrentLineFields();
	}

	public void close() throws AnalysisException {
		try {
			r.close();
		} catch (IOException e) {
			throw new AnalysisException("Error closing table file " + inTextStore.getPath() + ": "+e);				
		}			
	}
}
