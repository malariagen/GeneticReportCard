package org.cggh.common.fileIO;

import java.io.File;

import org.cggh.common.exceptions.*;
import org.cggh.common.textStore.*;

public class TableFileUtilities {
	
	public static void copyTableFile (File inFile, File outFile) throws AnalysisException {
		
		// Initialize output file
		OutputTextStore outStore =  isGzipped(outFile) ? new GzippedOutputTextStore(outFile.getParentFile(), getUngzippedFilename (outFile)) : new UncompressedOutputTextStore(outFile);
		BufferedTextOutput tableOut = new BufferedTextOutput(outStore, 1*1024*1024);
		StringBuffer sb = tableOut.getStringBuffer();
		
		// Initialize input files
		File inf = isGzipped(inFile) ? new File(inFile.getParentFile(), getUngzippedFilename (inFile)) : inFile;
		TableInput tableIn = new TableInput (new InputTextStore (inf));
		try {
			String[] fields = tableIn.getNextValidLine();
			sb.append(tableIn.getCurrentLine());
			while (true) {
				fields = tableIn.getNextValidLine();
				if (fields == null) {
					break;
				}
				sb.append('\n');
				sb.append(tableIn.getCurrentLine());
				tableOut.commitIfBufferFull();
			}
		} finally {
			tableIn.close();
			tableOut.close();
		}
	}
	

	public static void mergeTableFiles (File inFile1, File inFile2, File outFile) throws AnalysisException {
		mergeTableFiles (inFile1, inFile2, outFile, 1, null);
	}

	public static void mergeTableFiles (File inFile1, File inFile2, File outFile, int headerColumnCount) throws AnalysisException {
		mergeTableFiles (inFile1, inFile2, outFile, headerColumnCount, null);
	}

	public static void mergeTableFiles (File inFile1, File inFile2, File outFile, int headerColumnCount, String checkHeader) throws AnalysisException {

		// Initialize output file
		OutputTextStore outStore =  isGzipped(outFile) ? new GzippedOutputTextStore(outFile.getParentFile(), getUngzippedFilename (outFile)) : new UncompressedOutputTextStore(outFile);
		BufferedTextOutput tableOut = new BufferedTextOutput(outStore, 1*1024*1024);
		StringBuffer sb = tableOut.getStringBuffer();
		
		// Initialize input files
		File inf1 = isGzipped(inFile1) ? new File(inFile1.getParentFile(), getUngzippedFilename (inFile1)) : inFile1;
		File inf2 = isGzipped(inFile2) ? new File(inFile2.getParentFile(), getUngzippedFilename (inFile2)) : inFile2;
		TableInput tableInput1 = new TableInput (new InputTextStore (inf1));
		TableInput tableInput2 = new TableInput (new InputTextStore (inf2));
		
		int checkIdx = (checkHeader == null) ? -1 : tableInput1.getFieldIndex(checkHeader);
		
		// Write out the table header row
		mergeRowCells (tableInput1, tableInput2, sb, headerColumnCount);
		
		try {
			while (true) {
				String[] fields1 = tableInput1.getNextValidLine();
				String[] fields2 = tableInput2.getNextValidLine();
				if (fields1 == null) {
					break;
				}
				
				// Check the header columns if applicable
				if ((checkIdx >= 0) && !fields1[checkIdx].equals(fields2[checkIdx])) {
					throw new AnalysisException("At line number " + tableInput1.getLineNumber() + " there is a mismatch for header " + checkHeader 
							+ " ([" + fields1[checkIdx] + "] instead of [" + fields2[checkIdx] + "]) when merging files " 
							+ inFile1.getAbsolutePath() + " and " + inFile2.getAbsolutePath());
				}
				
				// Write out the merged row
				sb.append('\n');
				mergeRowCells (tableInput1, tableInput2, sb, headerColumnCount);
				tableOut.commitIfBufferFull();
			}
		} finally {
			tableInput1.close();
			tableInput2.close();
			tableOut.close();
		}
	}
	
	private static void mergeRowCells (TableInput tableInput1, TableInput tableInput2, StringBuffer sb, int headerColumnCount) {
		sb.append(tableInput1.getCurrentLine());
		String[] fields2 = tableInput2.getCurrentLineFields();
		for (int fIdx = headerColumnCount; fIdx < fields2.length; fIdx++) {
			sb.append('\t');
			sb.append(fields2[fIdx]);
		}
	}
	
	private static String getUngzippedFilename (File f) {
		String fname = f.getName();
		if (fname.endsWith(TextStore.GZIP_EXTENSION)) {
			return fname.substring(0, fname.length()-TextStore.GZIP_EXTENSION.length());
		}
		return fname;
	}
	
	private static boolean isGzipped (File f) {
		return f.getName().endsWith(TextStore.GZIP_EXTENSION);
	}
	
}
