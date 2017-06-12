package org.cggh.common.fileIO;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.fileIO.TabFileReader;
import org.cggh.common.textStore.InputTextStore;

import java.io.*;

public class ColumnFileReader {
	
	private TabFileReader tabReader;
	private String[]  columnNames;
	private String[]  currValues;
	
	public ColumnFileReader (InputTextStore inputStore) throws AnalysisException {

		if ((inputStore.getFile() == null) || !inputStore.getFile().exists()) {
			throw new AnalysisException("File "+inputStore.getPath()+" not found");
		}
		
		try {
			this.tabReader = new TabFileReader(inputStore.getReader());
		} catch (IOException e) {
			throw new AnalysisException("Error creating the TabFileReader: "+ e);
		}
		
		// Process headers
		String[] fields = null;
		try {
			fields = tabReader.getNextValidLine();
		} catch (IOException e) {
			throw new AnalysisException ("Error reading line #" + getLineNumber());
		}

		columnNames = new String[fields.length];
		for (int colIdx = 0; colIdx < columnNames.length; colIdx++) {
			columnNames[colIdx] = fields[colIdx].trim();
		}
		
		this.currValues = null;
	}
	
	public boolean nextRecord () throws AnalysisException {
		// Process fields
		try {
			currValues = tabReader.getNextValidLine();
		} catch (IOException e) {
			throw new AnalysisException ("Error reading line #" + getLineNumber());
		}
		
		if (currValues == null) {
			currValues = null;
			return false;
		}
		return true;
	}
	
	public String[] getValues () {
		return (currValues);
	}
	
	public String[] getColumnNames () {
		return (columnNames);
	}
	
	public int getLineNumber() {
		return tabReader.getLineNumber();
	}
	
	public void close () throws AnalysisException {
		if (tabReader != null) {
			try {
				tabReader.close();
			} catch (IOException e) {
				throw new AnalysisException("Error closing TabFileReader: "+ e);
			}
			
		}
	}
	
	public int[] getColumnIndexes (String[] queryColNames) throws AnalysisException {
		// Find indexes for the fields we will use
	    int[] queryColIdx = new int[queryColNames.length];
		for (int colIdx = 0; colIdx < queryColNames.length; colIdx++) {
			queryColIdx[colIdx] = getColumnIndex (queryColNames[colIdx]);
		}
		return queryColIdx;
	}
	
	public int getColumnIndex (String colName) throws AnalysisException {
		for (int fIdx = 0; fIdx < columnNames.length; fIdx++) {
			if (colName.equals(columnNames[fIdx])) {
				return fIdx;
			}
		}
		throw new AnalysisException("Column " + colName + " not found in file.");
	}
	
	public ColumnReader getColumnReader (String[] queryColNames) throws AnalysisException {
		return new ColumnReader (queryColNames);
	}

	public class ColumnReader {
	    protected String[]    queryColNames;
	    protected int[]       queryColIdx;

	    public ColumnReader (String[] queryColNames) throws AnalysisException {
		    this.queryColNames = queryColNames;
		    this.queryColIdx = getColumnIndexes (queryColNames);
	    }
	    
		public String[] getValues () {
			String[] queryValues = new String[queryColNames.length];
			for (int cIdx = 0; cIdx < queryColNames.length; cIdx++) {
				int vIdx = queryColIdx[cIdx];
				queryValues[cIdx] = currValues[vIdx];
			}
			return (queryValues);
		}
	}
}

