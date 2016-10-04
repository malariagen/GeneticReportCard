package org.cggh.bam;

import org.cggh.common.exceptions.*;
import org.cggh.common.threading.*;
import org.apache.commons.logging.*;
import java.io.*;


public class MultiSampleAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private Sample[] samples;
	private int      maxThreads = 0;
	
	public MultiSampleAnalysis (File sampleListFile, int maxThreads) throws AnalysisException  {
		this.samples = new SampleList(sampleListFile).getSamples();
		this.maxThreads = maxThreads;
	}
	
	public Sample[] getSamples() {
		return samples;
	}

	public void execute (SampleAnalysis analysis) throws AnalysisException, IOException  {

		// Analyze the samples one by one at all the loci		
		ParallelExecutableManager pem;
		if (maxThreads > 0) { 
			pem = new ParallelExecutableManager(maxThreads);
		} else {
			pem = new ParallelExecutableManager();
		}
		
		// Create one task per sample
		SampleAnalysisTask[] analysisTasks = new SampleAnalysisTask[samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			analysisTasks[sIdx] = new SampleAnalysisTask (analysis, samples[sIdx]);
			pem.addTask(analysisTasks[sIdx]);
		}

		// Analyze all the samples at all the loci
		pem.startExecution();
		pem.waitForCompletion();
		
		// Marge all the sample results files
		//mergeSampleResultFiles(samples);
	}
	
	private class SampleAnalysisTask implements Runnable {
		
		SampleAnalysis analysis;
		Sample sample;

		public SampleAnalysisTask (SampleAnalysis analysis, Sample sample) {
			this.analysis = analysis;
			this.sample = sample;
		}
		
		@Override
		public void run() {
			try {
				analysis.analyzeSample(sample);
			} catch (Exception e) {
				log.error("Error processing sample " + sample.getName() + ": "+e);
			}
		}
	}
}
