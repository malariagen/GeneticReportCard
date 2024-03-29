package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;


public abstract class BaseAnalysis {
	
	protected static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());

	protected File outRootFolder;
	
	public BaseAnalysis (File outRootFolder) throws AnalysisException  {
		// Create the root output directory
		this.outRootFolder = outRootFolder;
		if (!outRootFolder.exists()) {
			outRootFolder.mkdirs();
		}		
	}
	
	/* *************************************************************************
	 * Utility routined for merging the rows of many individual sample files
	 * *************************************************************************
	 */
	protected void mergeSampleDataFiles (File folder, String filenameSuffix, String[] fileHeaders, Sample[] samples, boolean warnIfMissing) throws AnalysisException {
		TableOutput out = new TableOutput (folder, "AllSamples_"+ filenameSuffix+".tab", fileHeaders, 64 * 1024);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleInFile = new File (folder, sample.getName()+"_"+filenameSuffix+".tab");
			appendTableFileContent (out, fileHeaders, sampleInFile, warnIfMissing);
		}
		out.close();
	}

	protected void appendTableFileContent (TableOutput mergedOut, String[] headers, File inFile, boolean warnIfMissing) throws AnalysisException {
		if (!inFile.exists()) {
			if (warnIfMissing) {
				log.warn("Sample file "+inFile.getAbsolutePath()+" not found.");
			}
			return;
		} else {			
			TableInput tif = new TableInput (new InputTextStore (inFile));
			// Check header consistency
			String[] fieldNames = tif.getFieldNames();
			if (fieldNames.length != (headers.length+1)) {
				log.error("Sample file "+inFile.getAbsolutePath()+": inconsistent headers - skipping.");
				tif.close();
				return;
			}
			for (int i = 0; i < headers.length; i++) {
				if (!fieldNames[i+1].equals(headers[i])) {
					log.error("Sample file "+inFile.getAbsolutePath()+": inconsistent headers - skipping.");
					tif.close();
					return;
				}
			}
			// Copy content
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					mergedOut.newRow();
					for (int i = 0; i < headers.length; i++) {
						mergedOut.appendValue(inFields[i+1]);
					}
				}
			} finally {
				tif.close();
			}
		}
	}
	
	
	/*
	 * Divide the samples into subfolders so we don't end up with thousands of files in the same folder
	 */
	protected File getSampleSubfolder (File rootFolder, Sample sample, boolean createIfMissing) {
		String batchName = sample.getBatch();
		String subFolderName;
		if (batchName == null) {
			subFolderName = sample.getName().substring(0, 4);
		} else {
			subFolderName = Sample.NO_BATCH.equals(batchName) ? "NO_BATCH" : batchName;
		}

		File subFolder;
		try {
			subFolder = FileUtilities.checkFolder(rootFolder, subFolderName, createIfMissing);
		} catch (AnalysisException e) {
			return null;
		}
		return subFolder;
	}
}
