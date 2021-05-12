package org.cggh.bam;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.fileIO.ColumnFileReader;
import org.cggh.common.textStore.InputTextStore;
import java.io.*;
import java.util.*;


public class SampleList {
	
	private Sample[] samples;
	private String[] batches;
	private File     sampleListFile;
	private boolean  checkBams;
	
	public SampleList (File sampleListFile) throws AnalysisException {
		this(sampleListFile, true);
	}
	
	public SampleList (File sampleListFile, boolean checkBams) throws AnalysisException {
		this.sampleListFile = sampleListFile;
		this.checkBams = checkBams;
		this.samples = readSamples();
		this.batches = extractBatches(samples);
	}
	
	public Sample[] getSamples() throws AnalysisException {
		return samples;
	}
	
	public String[] getBatches() throws AnalysisException {
		return batches;
	}
	
	private Sample[] readSamples() throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));
		String[] colNames = new String[]{"Batch","Sample","BamFile"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample(values[0], values[1], new File(values[2]));
			if (checkBams && !s.getBamFile().canRead()) {
				throw new AnalysisException("Error getting samples list: BAM file " + values[1] + " cannot be read.");
			}
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}

	private String[] extractBatches(Sample[] samples) throws AnalysisException {
		HashSet<String> bList = new HashSet<String>();
		for (Sample s : samples) {
			bList.add(s.getBatch());
		}
		String[] batches = bList.toArray(new String[bList.size()]);
		Arrays.sort(batches);
		return batches;
	}
}
