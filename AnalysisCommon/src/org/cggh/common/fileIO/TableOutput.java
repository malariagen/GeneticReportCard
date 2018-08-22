package org.cggh.common.fileIO;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;

import java.io.File;
import java.text.NumberFormat;


public class TableOutput {
		
	public static final int COMPRESSION_NONE = 1;
	public static final int COMPRESSION_GZIP = 2;
	
	public static final int LINE_NUMBERS_OFF = 1;
	public static final int LINE_NUMBERS_ON = 2;
	
	public static final int DEFAULT_MAX_DECIMAL_DIGITS = 4;

	protected static int defaultMaxDecimalDigits = DEFAULT_MAX_DECIMAL_DIGITS;
	public static void setDefaultMaxDecimalDigits (int maxFractionDigits) throws AnalysisException {
		defaultMaxDecimalDigits = maxFractionDigits;
	}
	
	protected boolean         initialized;
	protected boolean         lineNumbersOn;
	protected OutputTextStore outputStore;
	protected StringBuffer    sb;

	protected String[] colNames;
	protected int      rowCount;
	protected NumberFormat floatFormat;
	
	protected boolean firstFieldInFile;
	protected boolean firstFieldInLine;
	
	
	public TableOutput (File folder, String filename, String[] colNames, int bufferSize) throws AnalysisException {
		this (getOutputTextStore (folder, filename, COMPRESSION_NONE), colNames, bufferSize, LINE_NUMBERS_ON);
	}
	
	public TableOutput (File folder, String filename, String[] colNames, int bufferSize, int compression) throws AnalysisException {
		this (getOutputTextStore (folder, filename, compression), colNames, bufferSize, LINE_NUMBERS_ON);
	}
	
	public TableOutput (File folder, String filename, String[] colNames, int bufferSize, int compression, int lineNumbersOn) throws AnalysisException {
		this (getOutputTextStore (folder, filename, compression), colNames, bufferSize, lineNumbersOn);
	}
	
	public TableOutput (OutputTextStore outputStore, String[] colNames, int bufferSize) throws AnalysisException {
		this (outputStore, colNames, bufferSize, LINE_NUMBERS_ON);
	}
	
	public TableOutput (OutputTextStore outputStore, String[] colNames, int bufferSize, int lineNumbersOn) throws AnalysisException {
		this.initialized = false;
		this.lineNumbersOn = (lineNumbersOn == LINE_NUMBERS_ON);
		this.outputStore = outputStore;
		this.colNames    = colNames;
		this.sb          = new StringBuffer(bufferSize);
		
		setMaximumFractionDigits (defaultMaxDecimalDigits);
	}

	public void newRow () throws AnalysisException {
		newRow (null); // This method only works for the default TableOutput class, for which there is no row object (e.g. a SNP or sample)
	}
	
	public void newRow (Object rowObject) throws AnalysisException {
		
		// Do the headers if it's the first row
		if (!initialized) {
			FileUtilities.initializeFile(outputStore);
			firstFieldInFile = firstFieldInLine = true;
			if (lineNumbersOn) {
				appendFieldSeparator(sb);
				sb.append("Num");
			}
			writeRowHeaderColNames(sb);
			for (int colIdx = 0; colIdx < colNames.length; colIdx++) {
				appendFieldSeparator(sb);
				sb.append(colNames[colIdx]);
			}
			initialized = true;
		}
		
		// Start a new row, with the appropriate row headers
		commitIfBufferFull();
		rowCount++;
		firstFieldInLine = true;
		if (lineNumbersOn) {
			appendFieldSeparator(sb);
			sb.append(rowCount);
		}
		writeRowHeaderColValues (sb, rowObject);
	}
	
	protected void appendFieldSeparator(StringBuffer sb) {
		if (firstFieldInLine) {
			if (firstFieldInFile) {
				firstFieldInFile = false;
			} else {
				sb.append('\n');
			}
			firstFieldInLine = false;
		} else {
			sb.append('\t');
		}
	}
	
	// To be overridden
	protected void writeRowHeaderColNames (StringBuffer sb) {}
		
	// To be overridden
	protected void writeRowHeaderColValues (StringBuffer sb, Object rowObject) {}
		

	public void appendMultipleBlankValues (int colCount) throws AnalysisException {
		for (int i = 0; i < colCount; i++) {
			appendBlankValue ();
		}
	}
	
	public void appendBlankValue () throws AnalysisException {
		appendFieldSeparator(sb);
		sb.append('-');
	}
	
	public void appendValue (String value) throws AnalysisException {
		appendFieldSeparator(sb);
		if (value == null) {
			sb.append('-');
		} else {
			sb.append(value);
		}
	}
	
	public void appendMultipleValues (String[] values) throws AnalysisException {
		for (int i = 0; i < values.length; i++) {
			appendValue(values[i]);
		}
	}
	
	public void appendMultipleValues (int[] values) throws AnalysisException {
		for (int i = 0; i < values.length; i++) {
			appendValue(values[i]);
		}
	}
	
	public void appendValue (int value) throws AnalysisException {
		appendFieldSeparator(sb);
		if (value == Integer.MIN_VALUE) {
			sb.append('-');
		} else if (value == 0) {
			sb.append('0');
		} else {
			sb.append(value);
		}
	}
	
	public void appendValue (char value) throws AnalysisException {
		appendFieldSeparator(sb);
		if (value == Character.MIN_VALUE) {
			sb.append('-');
		} else {
			sb.append(value);
		}
	}
	
	public void appendValue (double value) throws AnalysisException {
		this.appendValue(value, floatFormat);
	}
	
	public void appendValue (double value, NumberFormat formatter) throws AnalysisException {
		appendFieldSeparator(sb);
		if (Double.isNaN(value)) {
			sb.append('-');
		} else if (value == 0.0) {
			sb.append('0');
		} else {
			sb.append(formatter.format(value));
		}
	}
	
	public void appendValue (boolean value) throws AnalysisException {
		appendFieldSeparator(sb);
		sb.append(Boolean.toString(value));
	}
	
	public void close() throws AnalysisException {
		commitIfHasContent();
		try {
			outputStore.closeWriter();
		} catch (Exception e) {
			throw new AnalysisException ("Error closing writer for "+outputStore.getPath()+":" + e);
		}
	}

	public void setFloatFormat (NumberFormat formatter) throws AnalysisException {
		this.floatFormat = formatter;
	}
		
	public void setMaximumFractionDigits (int maxFractionDigits) throws AnalysisException {
		floatFormat = NumberFormat.getNumberInstance();
		floatFormat.setMinimumFractionDigits(0);
		floatFormat.setMaximumFractionDigits(maxFractionDigits);
	}
	
	static private OutputTextStore getOutputTextStore (File folder, String filename, int compression) throws AnalysisException {
		OutputTextStore ots = null;
		switch(compression) {
		case COMPRESSION_NONE:
			ots = new UncompressedOutputTextStore(folder, filename);
			break;
		case COMPRESSION_GZIP:
			ots = new GzippedOutputTextStore(folder, filename);
			break;
		default:
			throw new AnalysisException("Error creating TableOutput: unknown compression flag: "+compression);
		}
		return ots;
	}
	
	public void commitIfBufferFull() throws AnalysisException {
		FileUtilities.commitIfBufferFull (sb, outputStore);	
	}

	public void commitIfHasContent() throws AnalysisException {
		FileUtilities.commitIfHasContent(sb, outputStore);	
	}
}
