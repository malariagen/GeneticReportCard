package org.cggh.bam;

import java.io.File;
import java.util.ArrayList;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.fileIO.ColumnFileReader;
import org.cggh.common.textStore.InputTextStore;

public class SampleList {
	
	private Sample[] samples;
	private File     sampleListFile;
	private boolean  checkBams;
	
	public SampleList (File sampleListFile) throws AnalysisException {
		this(sampleListFile, true);
	}
	
	public SampleList (File sampleListFile, boolean checkBams) throws AnalysisException {
		this.sampleListFile = sampleListFile;
		this.checkBams = checkBams;
		this.samples = readSamples();
	}
	
	public Sample[] getSamples() throws AnalysisException {
		return samples;
	}
	
	private Sample[] readSamples() throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));
		String[] colNames = new String[]{"Sample","BamFile"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample(values[0], new File(values[1]));
			if (checkBams && !s.getBamFile().canRead()) {
				throw new AnalysisException("Error getting samples list: BAM file " + values[1] + " cannot be read.");
			}
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}

}
