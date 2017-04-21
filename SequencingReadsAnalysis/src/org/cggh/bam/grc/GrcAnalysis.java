package org.cggh.bam.grc;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.bam.target.SampleAnalyzer.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class GrcAnalysis extends SampleTargetAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private GrcConfig config;
	
	public GrcAnalysis (File configFile, File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, chrMapFile, outRootFolder);
		
		// Parse configuration file
		config = new GrcConfig (configFile);
		registerLoci(config.getLoci());
	}

	/* **********************************************************************
	 * Single sample processing
	 * **********************************************************************
	 */
	public void analyzeSample (Sample sample) throws AnalysisException  {
		log.info("Starting " + sample.getName());
		try {
			SampleAnalyzer analyzer = new SampleAnalyzer (config, sample);
			SampleResults sr = analyzer.analyzeSample();
			
			// Write out the results
			outputSampleResults (sr);
			
		} catch (Exception e) {
			log.info("Aborting " + sample.getName());
			log.error("Error processing BAM file for sample "+ sample.getName() + ": "+e);
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}
	
	private static final String[] LOCUS_COUNTS_HEADERS = new String[] {"Sample","Locus","Aligned","Misaligned"};
	private static final String[] ALLELE_COUNTS_HEADERS = new String[] {"Sample","Locus","Target","Allele","Amino","Count"};
	
	public void outputSampleResults (SampleResults sr) throws AnalysisException, IOException  {
		Sample sample = sr.getSample();
		File outFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
		TableOutput locusOut = new TableOutput (outFolder, sample.getName()+".locusCoverage.tab", LOCUS_COUNTS_HEADERS, 64 * 1024);		
		TableOutput alleleOut = new TableOutput (outFolder, sample.getName()+".alleles.tab", ALLELE_COUNTS_HEADERS, 64 * 1024);		

		SampleLocusResult[] locusResults = sr.getLocusResults();
		for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
			SampleLocusResult locusResult = locusResults[lIdx];
			TargetLocus locus = locusResult.getLocus();
			locusOut.newRow();
			locusOut.appendValue(sample.getName());
			locusOut.appendValue(locus.getName());
			locusOut.appendValue(locusResult.getAlignedCount());
			locusOut.appendValue(locusResult.getMisalignedCount());
			
			SampleTargetResult[] targetResults = locusResult.getTargetResults();
			for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
				SampleTargetResult targetResult = targetResults[tIdx];
				Target target = targetResult.getTarget();
				
				LabelCounter[] alleleCounters = targetResult.getNtAlleleCounters().getSortedCounters();
				for (int aIdx = 0; aIdx < alleleCounters.length; aIdx++) {
					alleleOut.newRow();
					alleleOut.appendValue(sample.getName());
					alleleOut.appendValue(locus.getName());
					alleleOut.appendValue(target.getName());
					
					String alleleSeq = alleleCounters[aIdx].getLabel();
					alleleOut.appendValue(alleleSeq);
					alleleOut.appendValue(SequenceUtilities.translateNtSequence(alleleSeq));
					alleleOut.appendValue(alleleCounters[aIdx].getCount());
				}
			}
		}
		locusOut.close();
		alleleOut.close();
	}
	
	
	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	public void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {
		
		// Go through all the samples, reading in all the allele counts for each target
		SampleTargetResult[][] allTargetResults = readAllSampleResults (samples);
		
		// Organize the coverage info by locus from sample-wise coverage data files
		processLocusCoverageInfo (samples);

		// Analyze one target at a time
		SampleCaller caller = new SampleCaller(5, 2);
		AminoSampleCall[][] allCalls = new AminoSampleCall[allTargets.length][samples.length];			
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			Target target = allTargets[tIdx];
			String targetName = allTargetNames[tIdx];
			
			// Get results for all the samples at this target
			SampleTargetResult[] tSampleResults = allTargetResults[tIdx];
			
			// Count the samples per allele in a master lists
			LabelCounters alleleSampleCounters = new LabelCounters();
			for (int sIdx = 0; sIdx < tSampleResults.length; sIdx++) {
				// Remove singleton alleles for the sample
				SampleTargetResult sampleResult = tSampleResults[sIdx];
				if (sampleResult == null) {
					continue;
				}
				sampleResult.cleanupTargetAlleles();
				LabelCounter[] ac = sampleResult.getNtAlleleCounters().getSortedCounters();
				for (int aIdx = 0; aIdx < ac.length; aIdx++) {
					alleleSampleCounters.increment(ac[aIdx].getLabel());						
				}
			}

			// Write out the Sample Summary for this target
			outputSampleTargetSummary (tSampleResults, outRootFolder, "SummaryBySample."+targetName+".tab");

			// Compute some target Allele statistics
			AlleleStats[] alleleStats = computeTargetStats (tSampleResults, alleleSampleCounters);

			// Write out pre-filtered allele read count table and allele stats summary table
			outputAlleleReadCounts (tSampleResults, alleleStats, outRootFolder, "AlleleSampleCount.beforeMinAlt."+targetName+".tab");
			outputAlleleStats (alleleStats, outRootFolder, "AlleleStats.beforeMinAlt."+targetName+".tab");

			// Now go through and filter the alleles to reject likely artefact
			// The test is MinAlt: homozygous or >=10 reads in at least one sample
			ArrayList<AlleleStats> filteredAllelelStatList = new ArrayList<AlleleStats>();
			for (int aIdx = 0; aIdx < alleleStats.length; aIdx++) {
				if ((alleleStats[aIdx].maxReadsFraction == 1.0) || (alleleStats[aIdx].maxReads >= 10)) {
					filteredAllelelStatList.add(alleleStats[aIdx]);
				}
			}
			
			AlleleStats[] filteredAlleleStats = filteredAllelelStatList.toArray(new AlleleStats[filteredAllelelStatList.size()]);
			String[] filteredAlleles = new String[filteredAlleleStats.length];
			for (int aIdx = 0; aIdx < filteredAlleleStats.length; aIdx++) {
				filteredAlleles[aIdx] = filteredAlleleStats[aIdx].ntSequence;
			}
			
			// Go through the samples, keeping only the alleles that passed the MinAlt filter.
			LabelCounterFilter alleleListFilter = new LabelCounterFilterByLabelList (filteredAlleles);
			for (int sIdx = 0; sIdx < tSampleResults.length; sIdx++) {
				tSampleResults[sIdx].filterAlleleCounters(alleleListFilter);
			}

			// Write out an allele read count table and an allele stats summary table
			outputAlleleReadCounts (tSampleResults, filteredAlleleStats, outRootFolder, "AlleleSampleCount.final."+targetName+".tab");
			outputAlleleStats (filteredAlleleStats, outRootFolder, "AlleleStats.final."+targetName+".tab");
			
			// Final step: make a call 
			for (int sIdx = 0; sIdx < tSampleResults.length; sIdx++) {
				SampleCall ntCall = caller.callSample(tSampleResults[sIdx]);
				allCalls[tIdx][sIdx] = new AminoSampleCall(ntCall, target);
			}
			// Write out the final calls
			outputSampleCalls (allCalls[tIdx], samples, target, outRootFolder, "CallsBySample.final."+targetName+".tab");
			allTargets[tIdx] = target;
		}
		
		// Finally write out the overall results table
		outputSampleCallsAllTargets (allTargets, samples, allCalls);
	}


	
	private SampleTargetResult[][] readAllSampleResults (Sample[] samples) throws AnalysisException {

		SampleTargetResult[][] targetResults = new SampleTargetResult[allTargets.length][samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
				targetResults[tIdx][sIdx] = new SampleTargetResult(allTargets[tIdx], samples[sIdx], new LabelCounters());
			}
		}
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
			File sampleFile = new File (sampleFolder, sample.getName()+".alleles.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFIdx       = tif.getFieldIndex("Locus");
			int targetFIdx      = tif.getFieldIndex("Target");
			int alleleFldIdx    = tif.getFieldIndex("Allele");
			int countFldIdx     = tif.getFieldIndex("Count");
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					//String sampleName = inFields[sampleFIdx];
					String tName = inFields[locusFIdx]+"_"+inFields[targetFIdx];
					int tIdx = getTargetIndex (tName);
					
					String allele = inFields[alleleFldIdx];
					int count = Integer.parseInt(inFields[countFldIdx]);
					targetResults[tIdx][sIdx].getNtAlleleCounters().setCount(allele, count);
				}
			} finally {
				tif.close();
			}
		}
		return targetResults;
	}
	
	public void processLocusCoverageInfo (Sample[] samples) throws AnalysisException, IOException  {
		
		int[][] alignedCount = new int[loci.length][samples.length];
		int[][] misalignedCount = new int[loci.length][samples.length];
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
			File sampleFile = new File (sampleFolder, sample.getName()+".locusCoverage.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFldIdx      = tif.getFieldIndex("Locus");
			int alignedFldIdx    = tif.getFieldIndex("Aligned");
			int misalignedFldIdx = tif.getFieldIndex("Misaligned");
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					//String sampleName = inFields[sampleFIdx];
					String locusName = inFields[locusFldIdx];
					int lIdx = getLocusIndex (locusName);
					alignedCount[lIdx][sIdx] = Integer.parseInt(inFields[alignedFldIdx]);
					misalignedCount[lIdx][sIdx] = Integer.parseInt(inFields[misalignedFldIdx]);
				}
			} finally {
				tif.close();
			}
		}
		
		// Write out the results by locus
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus locus = loci[lIdx];
			String locusName = locusNames[lIdx];
			TableOutput locusOut = new TableOutput (outRootFolder, "LocusCoverage."+locusName+".tab", LOCUS_COUNTS_HEADERS, 64 * 1024);
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				locusOut.newRow();
				locusOut.appendValue(samples[sIdx].getName());
				locusOut.appendValue(locus.getName());
				locusOut.appendValue(alignedCount[lIdx][sIdx]);
				locusOut.appendValue(misalignedCount[lIdx][sIdx]);
			}
			locusOut.close();
		}
	}
	
	private void outputSampleCalls (AminoSampleCall[] targetCalls, Sample[] samples, Target target, File outFolder, String filename) throws AnalysisException {
		int bufferSize = 64 * 1024;
		String[] headers = new String[] {"Sample","Call","Alleles","AlleleReads","NtAlleles","NtAlleleReads"};
		TableOutput sampleCallsOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < targetCalls.length; sIdx++) {
			AminoSampleCall call = targetCalls[sIdx];
			boolean isLenient = call.isLenient();
			String callString = isLenient ? "["+call.getCallString()+"]" : call.getCallString();
			sampleCallsOut.newRow();
			sampleCallsOut.appendValue(samples[sIdx].getName());
			sampleCallsOut.appendValue(callString);
			sampleCallsOut.appendValue(call.getAminoAllele());
			sampleCallsOut.appendValue(call.getAminoAlleleSummary());
			sampleCallsOut.appendValue(call.getAllele());
			sampleCallsOut.appendValue(call.getAlleleSummary());
		}
		sampleCallsOut.close();
	}
	

    private void outputSampleCallsAllTargets (Target[] allTargets, Sample[] samples, AminoSampleCall[][] allCalls) throws AnalysisException {

    	int bufferSize = 1024 * 1024;
    	// Get the file headers
		ArrayList<String> headerList = new ArrayList<String>();
		headerList.add("Sample");
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			Target target = allTargets[tIdx];
			String targetAminoRefSeq = SequenceUtilities.translateNtSequence(target.getTargetRefSeq());
			headerList.add(target.getName()+"["+targetAminoRefSeq+"]");
		}
		String[] headers = headerList.toArray(new String[headerList.size()]);

		TableOutput sampleCallsOut = new TableOutput (outRootFolder, "AllCallsBySample.tab", headers, bufferSize);
		TableOutput sampleNrefCallsOut = new TableOutput (outRootFolder, "AllCallsNrefBySample.tab", headers, bufferSize);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			sampleCallsOut.newRow();
			sampleCallsOut.appendValue(samples[sIdx].getName());
			
			sampleNrefCallsOut.newRow();
			sampleNrefCallsOut.appendValue(samples[sIdx].getName());
			
			for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
				AminoSampleCall call = allCalls[tIdx][sIdx];
				String aaAllele = call.getAminoAllele();
				String aaNrefAllele = call.getAminoNrefAllele();
				if (!call.isMissing() && call.isLenient()) {
					aaAllele = aaAllele.toLowerCase();
					aaNrefAllele = aaNrefAllele.toLowerCase();							
				}
				sampleCallsOut.appendValue(aaAllele);
				sampleNrefCallsOut.appendValue(aaNrefAllele);
			}
		}
		sampleCallsOut.close();
		sampleNrefCallsOut.close();
	}

	
	/* ==========================================================
	 * Allele statistics computation
	 * ==========================================================
	 */
	public class AlleleStats {
		String ntSequence;
		String aaSequence;
		String displayLabel;
		int    sampleCount;
		int    maxReads;
		double maxReadsFraction;
		
		public AlleleStats(String ntSequence) {
			super();
			this.ntSequence = ntSequence;
			this.aaSequence = SequenceUtilities.translateNtSequence(ntSequence);
			this.displayLabel = ntSequence+"["+aaSequence+"]";
		}
	}	
	
	private AlleleStats[] computeTargetStats (SampleTargetResult[] sampleResults, LabelCounters alleleSampleCounters) throws AnalysisException {
		LabelCounter[] alleleCounters = alleleSampleCounters.getSortedCounters();
		AlleleStats[] aStats = new AlleleStats[alleleCounters.length];
		for (int aIdx = 0; aIdx < alleleCounters.length; aIdx++) {
			aStats[aIdx] = new AlleleStats(alleleCounters[aIdx].getLabel());
		}
		
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			LabelCounters aCounters = sampleResults[sIdx].getNtAlleleCounters();
			double totalSampleReads = aCounters.getTotal();
			for (int aIdx = 0; aIdx < alleleCounters.length; aIdx++) {
				AlleleStats stat = aStats[aIdx];
				LabelCounter c = aCounters.getCounter(stat.ntSequence);
				if (c != null) {
					int readCount = c.getCount();
					double readsFraction = ((double)readCount) / totalSampleReads;
					stat.sampleCount++;
					stat.maxReads = (readCount > stat.maxReads) ? readCount : stat.maxReads;
					stat.maxReadsFraction = (readsFraction > stat.maxReadsFraction) ? readsFraction : stat.maxReadsFraction;
				}
			}
		}
		return aStats;
	}

	private void outputAlleleReadCounts (SampleTargetResult[] sampleResults, AlleleStats[] aStats, File outFolder, String filename) throws AnalysisException {
		int bufferSize = 64 * 1024;
		String[] headers = new String[aStats.length+1];
		headers[0] = "Sample";
		for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
			headers[aIdx+1] = aStats[aIdx].displayLabel;
		}
		TableOutput alleleReadsOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			SampleTargetResult sampleResult = sampleResults[sIdx];
			alleleReadsOut.newRow();
			alleleReadsOut.appendValue(sampleResult.getSample().getName());
			for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
				LabelCounter c = sampleResult.getNtAlleleCounters().getCounter(aStats[aIdx].ntSequence);
				alleleReadsOut.appendValue((c != null) ? c.getCount() : 0);
			}
		}
		alleleReadsOut.close();
	}

	private void outputAlleleStats (AlleleStats[] hStats, File outFolder, String filename) throws AnalysisException {
		int bufferSize = 128 * 1024;
		String[] headers = new String[] {"Allele","SampleCount","MaxReads","MaxReadFraction"};
		TableOutput alleleStatsOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			alleleStatsOut.newRow();
			alleleStatsOut.appendValue(hStats[hIdx].displayLabel);
			alleleStatsOut.appendValue(hStats[hIdx].sampleCount);
			alleleStatsOut.appendValue(hStats[hIdx].maxReads);
			alleleStatsOut.appendValue(hStats[hIdx].maxReadsFraction);
		}
		alleleStatsOut.close();
	}

	private void outputSampleTargetSummary (SampleTargetResult[] sampleResults, File outFolder, String filename) throws AnalysisException {
		// Write out the Sample Summary for this target
		int bufferSize = 16 * 1024;
		String[] headers = new String[]{"Sample","Alleles","AlleleCount"};
		TableOutput targetSummaryOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			SampleTargetResult sResult = sampleResults[sIdx];
			targetSummaryOut.newRow();
			targetSummaryOut.appendValue(sResult.getSample().getName());
			targetSummaryOut.appendValue(sResult.getNtAlleleCounters().getSummary());
			targetSummaryOut.appendValue(sResult.getNtAlleleCounters().getSize());
		}
		targetSummaryOut.close();
	}


	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		
		public static void main(String[] args) {
			if (args.length < 7) {
				log.error("Usage: org.cggh.bam.grc.GrcAnalysis$SingleSample <configFile> <sampleName> <bamFile> <chrMap> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			String sampleId = args[1];					log.info("SampleId: "+sampleId);
			File sampleBamFile = new File(args[2]);	    log.info("SampleBamFile: "+sampleBamFile.getAbsolutePath());
			String sampleChrMapName = args[3];			log.info("SampleChrMapName: "+sampleChrMapName);
			File refFastaFile = new File(args[4]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[5]);		log.info("ChrMapFile: "+chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[6]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample sample = new Sample (sampleId, sampleBamFile, sampleChrMapName);
				GrcAnalysis task = new GrcAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
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
			if (args.length < 5) {
				log.error("Usage: org.cggh.bam.grc.GrcAnalysis$MultiSample <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);		log.info("ChrMapFile: "+chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			
			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			
			try {	
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				GrcAnalysis task = new GrcAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
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
			if (args.length < 5) {
				log.error("Usage: org.cggh.bam.grc.GrcAnalysis$MergeResults <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);		log.info("ChrMapFile: "+chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				GrcAnalysis task = new GrcAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
