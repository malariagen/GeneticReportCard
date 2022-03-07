package org.cggh.bam;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.fileIO.ColumnFileReader;
import org.cggh.common.textStore.InputTextStore;
import java.io.*;
import java.util.*;


public class SampleList {
	
	private Sample[] samples;
	private File     sampleListFile;
	
	public SampleList (File sampleListFile) throws AnalysisException {
		this(sampleListFile, true);
	}
	
	public SampleList (File sampleListFile, boolean checkBams) throws AnalysisException {
		this.sampleListFile = sampleListFile;
		this.samples = readSamples(checkBams);
	}
	
	public Sample[] getSamples() throws AnalysisException {
		return samples;
	}
	
	private Sample[] readSamples(boolean  checkBams) throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));

		boolean hasBatch = cfr.hasColumn("Batch");
		String[] colNames = hasBatch ? new String[]{"Batch","Sample","BamFile"} : new String[]{"Sample","BamFile"};
		
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();

			String batchName = null;
			if (hasBatch) {
				batchName = values[0].trim();
				if (batchName.isEmpty()) {
					batchName = Sample.NO_BATCH;
				}
			}
			String sampleName = hasBatch ? values[1] : values[0];
			String bamFilename = hasBatch ? values[2] : values[1];
			
			Sample s = new Sample(batchName, sampleName, new File(bamFilename));
			if (checkBams && !s.getBamFile().canRead()) {
				throw new AnalysisException("Error getting samples list: BAM file " + values[1] + " cannot be read.");
			}
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}
}
