package org.cggh.common.fileIO;

import org.cggh.common.exceptions.*;
import org.cggh.common.textStore.*;
import java.io.*;


public class SelectableColumnTableInput extends TableInput {
	
	private String[] selColumnNames;
	private int[]    selColumnIndexes;

	public SelectableColumnTableInput (File inFile, String[] columnNames) throws AnalysisException {
		this(new InputTextStore (inFile.getParentFile(), inFile.getName()), columnNames);
	}
	
	public SelectableColumnTableInput (InputTextStore inTextStore, String[] columnNames) throws AnalysisException {
		this (inTextStore, "\t", DelimitedReader.DEFAULT_COMMENT_PREFIX, columnNames);
	}
	
	public SelectableColumnTableInput (File inFile, String delimiter, String commentPrefix, String[] columnNames) throws AnalysisException {
		this(new InputTextStore (inFile.getParentFile(), inFile.getName()), delimiter, commentPrefix, columnNames);
	}
	
	public SelectableColumnTableInput (InputTextStore inTextStore, String delimiter, String commentPrefix, String[] columnNames) throws AnalysisException {
		super(inTextStore, delimiter, commentPrefix);
		this.selColumnNames = columnNames;
		this.selColumnIndexes = new int[columnNames.length];
		
		for (int i = 0; i < selColumnNames.length; i++) {
			String cName = selColumnNames[i];
			int idx = -1;
			for (int fIdx = 0; fIdx < fieldNames.length; fIdx++) {
				if (fieldNames[fIdx].equals(cName)) {
					idx = fIdx;
					break;
				}
			}
			if (idx < 0) {
				throw new AnalysisException("Error opening table file "+ inTextStore.getPath()+ ": column '"+ cName+ "' not found");
			}
			selColumnIndexes[i] = idx;
		}
	}
	
	@Override
	public String[] getFieldNames () {
		return selColumnNames;
	}
	
	@Override
	public int getFieldIndex (String fieldName) {
		for (int i = 0; i < selColumnNames.length; i++) {
			if (fieldName.equals(selColumnNames[i])) {
				return i;
			}
		}
		return -1;
	}
	
	@Override
	public String[] getNextValidLine() throws AnalysisException {
		try {
			String[] values = r.getNextValidLine();
			if (values == null) {
				return null;
			}
			String[] selValues = new String[selColumnNames.length];
			for (int i = 0; i < selColumnNames.length; i++) {
				int selColumnIdx = selColumnIndexes[i];
				selValues[i] = values[selColumnIdx];
			}
			return selValues;
			
		} catch (IOException e) {
			throw new AnalysisException("Error reading record from file " + inTextStore.getPath() + " at line " + r.getLineNumber() + ": "+e);				
		}			
	}
}
