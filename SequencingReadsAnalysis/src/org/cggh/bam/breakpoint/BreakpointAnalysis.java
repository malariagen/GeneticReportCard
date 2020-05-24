package org.cggh.bam.breakpoint;

import org.cggh.bam.*;
import org.cggh.bam.breakpoint.SampleBreakpointAnalyzer.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class BreakpointAnalysis extends SampleAnalysis  {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private BreakpointConfig             config;
	
	
	public BreakpointAnalysis (File configFile, File outRootFolder) throws AnalysisException  {
		super (outRootFolder);

		// Parse configuration file
		config = new BreakpointConfig (configFile);
	}

	/* **********************************************************************
	 * Single sample processing
	 * **********************************************************************
	 */
	public void analyzeSample (Sample sample) throws AnalysisException  {
		log.info("Starting " + sample.getName());
		try {
			File outFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
			SampleBreakpointAnalyzer analyzer = new SampleBreakpointAnalyzer (config, sample, outFolder);
			SampleResult sr = analyzer.analyzeSample();
			
			// Analyze unlisted alleles: if they are very similar to listed ones, and can be 
			// assigned to a sample class, add them to the relevant counter
			//analyzeUnlistedAlleles (sr);
			
			// Write out the results
			outputSampleResults (sr, outFolder);
			
		} catch (Exception e) {
			log.info("Aborting " + sample.getName());
			log.error("Error processing BAM file for sample "+ sample.getName() + ": "+e);
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}
	
	//private static final String[] LISTED_ALLELES_HEADERS = {"Sample","Locus","Target","Allele","Count","SampleClass"};
	//private static final String[] UNLISTED_ALLELES_HEADERS = {"Sample","Locus","Target","Allele","Count","Proportion","Closest","Diff"};
	private static final String[] FLANK_HEADERS = {"Sample","Breakpoint","RemainderSeq","MatchSeq","FlankSeq"};
	
	
	/*
	 * Write out the results for this sample into two files: one of counts of listed sample class-specific alleles,
	 * and one for alleles that were not listed.
	 */
	protected void outputSampleResults (SampleResult sr, File outFolder) throws AnalysisException, IOException  {
		Sample sample = sr.sample;

		TableOutput flankOut = new TableOutput (outFolder, sample.getName()+".rbpFlanks.tab", FLANK_HEADERS, 64 * 1024);		
		ReadResult[] rResults = sr.rResults;
		for (int rrIdx = 0; rrIdx < rResults.length; rrIdx++) {
			ReadResult rr = rResults[rrIdx];
			flankOut.newRow();
			flankOut.appendValue(sample.getName());
			flankOut.appendValue(rr.breakpoint.id);
			flankOut.appendValue(rr.remainderSeq);
			flankOut.appendValue(rr.matchSeq);
			flankOut.appendValue(rr.flankSeq);
		}
		flankOut.close();
	}
	

	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	protected void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {
		mergeResultFiles(samples, "rbpFlanks");
	}

	private void mergeResultFiles(Sample[] samples, String filenameSuffix) throws AnalysisException, IOException {
		//TableOutput mergeOut = new TableOutput(outRootFolder, "AllSamples." + locus.getName() + filenameSuffix + ".tab", fieldHeaders, 1048576);
		TableOutput mergeOut = null;
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
			File sampleFile = new File(sampleFolder, String.valueOf(sample.getName()) + "." + filenameSuffix + ".tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
			} else {
				TableInput tif = new TableInput(sampleFile);
				// Get the headers, removing the "Num" automatic field at the start
				String[] fieldHeaders = tif.getFieldNames();
				fieldHeaders = Arrays.copyOfRange(fieldHeaders, 1, fieldHeaders.length); 
				if (mergeOut == null) {
					mergeOut = new TableOutput(outRootFolder, "AllSamples." + filenameSuffix + ".tab", fieldHeaders, 64 * 1024);
				}
				try {
					String[] inFields;
					while ((inFields = tif.getNextValidLine()) != null) {
						mergeOut.newRow();
						for (int fIdx = 0; fIdx < fieldHeaders.length; fIdx++) {
							mergeOut.appendValue(inFields[fIdx + 1]);
						}
					}
				} finally {
					tif.close();
				}
			}
		}
		if (mergeOut != null) {
			mergeOut.close();
		}
	}

	
	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		
		public static void main(String[] args) {
			if (args.length < 4) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$SingleSample <configFile> <sampleName> <bamFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			String sampleId = args[1];					log.info("SampleId: "+sampleId);
			File sampleBamFile = new File(args[2]);	    log.info("SampleBamFile: "+sampleBamFile.getAbsolutePath());
			File rootFolder = new File(args[3]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				BreakpointAnalysis task = new BreakpointAnalysis(configFile, rootFolder);
				Sample sample = new Sample (sampleId, sampleBamFile);
				task.analyzeSample(sample);	
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
	
	/* ==========================================================
	 * Multi Sample Execution
	 * ==========================================================
	 */
	public static class MultiSample {
		public static void main(String[] args) {
			if (args.length < 3) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MultiSample <configFile> <sampleListFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File rootFolder = new File(args[2]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			
			try {	
				BreakpointAnalysis task = new BreakpointAnalysis(configFile, rootFolder);
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				multi.execute(task);
				task.analyzeAllSampleResults(multi.getSamples());
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
	
	/* ==========================================================
	 * Sample Result Merging
	 * ==========================================================
	 */
	public static class MergeResults {
		public static void main(String[] args) {
			if (args.length < 3) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MergeResults <configFile> <sampleListFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File rootFolder = new File(args[2]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				BreakpointAnalysis task = new BreakpointAnalysis(configFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
