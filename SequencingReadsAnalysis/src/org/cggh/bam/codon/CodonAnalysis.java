package org.cggh.bam.codon;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.bam.target.SampleAnalyzer.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.cggh.common.util.FileUtilities;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class CodonAnalysis extends SampleTargetAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	public final static boolean APPLY_MIN_ALT = false;
	
	private CodonConfig config;
	
	public CodonAnalysis (File configFile, File refFastaFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, outRootFolder);
		
		// Parse configuration file
		config = new CodonConfig (configFile);
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
			String sampleName = sample.getName();
			String excMsg = e.toString();
			log.info("Aborting " + sampleName);
			log.error("Error processing BAM file for sample "+ sampleName + ": "+excMsg);
			try {
				recordSampleFailure (sampleName, excMsg);
			} catch (IOException e1) {}
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}
	
	private synchronized void recordSampleFailure (String sampleName, String excMsg) throws AnalysisException, IOException {
		File errorFile = new File (outRootFolder, "FailedSamples.tab");
		if (!errorFile.exists()) {
			FileUtilities.writeFileContent("Sample\tError", errorFile);
		}
		FileUtilities.appendFileContent("\n"+sampleName+"\t"+excMsg, errorFile);
	}
	
	private static final String[] LOCUS_COUNTS_HEADERS = new String[] {"Batch","Sample","Locus","Target","Aligned","Misaligned","Covering","Calls","LowQuality"};
	private static final String[] CALL_HEADERS = new String[] {"Batch","Sample","Locus","Target","Call","Amino","AminoNref","Nt","NtNref","Counts"};
	private static final String[] ALLELE_COUNTS_HEADERS = new String[] {"Batch","Sample","Locus","Target","Allele","Amino","Count"};
	
	public void outputSampleResults (SampleResults sr) throws AnalysisException, IOException  {
		Sample sample = sr.getSample();
		File outFolder = getSampleSubfolder (outRootFolder, sample, true);
		TableOutput locusOut = new TableOutput (outFolder, sample.getName()+".locusCoverage.tab", LOCUS_COUNTS_HEADERS, 64 * 1024);
		TableOutput callsOut = new TableOutput (outFolder, sample.getName()+".calls.tab", CALL_HEADERS, 64 * 1024);
		TableOutput alleleOut = new TableOutput (outFolder, sample.getName()+".alleles.tab", ALLELE_COUNTS_HEADERS, 64 * 1024);

		SampleLocusResult[] locusResults = sr.getLocusResults();
		for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
			SampleLocusResult locusResult = locusResults[lIdx];
			TargetLocus locus = locusResult.getLocus();
			SampleTargetResult[] targetResults = locusResult.getTargetResults();
			for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
				SampleTargetResult targetResult = targetResults[tIdx];
				Target target = targetResult.getTarget();
				
				// Output target coverage data (numbers of reads retrieves, used, discarded, etc.)
				locusOut.newRow();
			    locusOut.appendValue(sample.getBatch());
			    locusOut.appendValue(sample.getName());
			    locusOut.appendValue(locus.getName());
			    locusOut.appendValue(target.getName());
				locusOut.appendValue(locusResult.getAlignedCount());
			    locusOut.appendValue(locusResult.getMisalignedCount());
			    int hqCount = targetResult.getNtAlleleCounters().getTotal();
			    int lqCount = targetResult.getLowQualityCount();
			    locusOut.appendValue(hqCount+lqCount);
			    locusOut.appendValue(hqCount);
			    locusOut.appendValue(lqCount);
				
				// Output target calls (nt and amino)
			    callsOut.newRow();
			    callsOut.appendValue(sample.getBatch());
			    callsOut.appendValue(sample.getName());
			    callsOut.appendValue(locus.getName());
			    callsOut.appendValue(target.getName());
			    
			    callsOut.appendValue(targetResult.getAminoCall().getCallString());
			    callsOut.appendValue(targetResult.getAminoCall().getAminoAllele());
				callsOut.appendValue(targetResult.getAminoCall().getAminoNrefAllele());
			    callsOut.appendValue(targetResult.getNtCall().getAllele());
				callsOut.appendValue(targetResult.getNtCall().getNrefAllele());
				callsOut.appendValue(targetResult.getAminoCall().getAminoAlleleSummary());
				
				// Output allele read count detailed data
				LabelCounter[] alleleCounters = targetResult.getNtAlleleCounters().getSortedCounters();
				for (int aIdx = 0; aIdx < alleleCounters.length; aIdx++) {
					alleleOut.newRow();
					alleleOut.appendValue(sample.getBatch());
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
		callsOut.close();
		alleleOut.close();
	}
	
	
	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	public void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {
		
		// Go through all the samples, reading in all the allele counts and calls for each target
		SampleTargetResult[][] allTargetResults = readAllSampleResults (samples);
		AminoSampleCall[][]    allTargetCalls = readAllSampleCalls (samples);
		
		// Organize the coverage info by locus from sample-wise coverage data files
		processLocusCoverageInfo (samples);

		// Analyze one target at a time
		//AminoSampleCall[][] allCalls = new AminoSampleCall[allTargets.length][samples.length];			
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
				// Remove sequences with undetermined nucleotide, and singleton reads			
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

			// Write out an allele read count table and an allele stats summary table
			outputAlleleReadCounts (tSampleResults, alleleStats, outRootFolder, "AlleleSampleCount."+targetName+".tab");
			outputAlleleStats (alleleStats, outRootFolder, "AlleleStats."+targetName+".tab");
		
			// Write out the final calls
			outputSampleCalls (allTargetCalls[tIdx], samples, target, outRootFolder, "CallsBySample."+targetName+".tab");
			allTargets[tIdx] = target;
		}
		
		// Finally write out the overall results table
		outputSampleCallsAllTargets (allTargets, samples, allTargetCalls);
	}

	
	private AminoSampleCall[][] readAllSampleCalls (Sample[] samples) throws AnalysisException {
		// Go through all the samples, reading in all the calls for each target, and put them in the amino call objects
		AminoSampleCall[][] targetCalls = new AminoSampleCall[allTargets.length][samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
			File sampleFile = new File (sampleFolder, sample.getName()+".calls.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
		    
			int locusFIdx    = tif.getFieldIndex("Locus");
			int targetFIdx   = tif.getFieldIndex("Target");
			int callFldIdx   = tif.getFieldIndex("Call");
			int aaFldIdx     = tif.getFieldIndex("Amino");
			int aaNrefFldIdx = tif.getFieldIndex("AminoNref");
			int ntFldIdx     = tif.getFieldIndex("Nt");
			int ntNrefFldIdx = tif.getFieldIndex("NtNref");
			int countsFldIdx = tif.getFieldIndex("Counts");

			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					//String sampleName = inFields[sampleFIdx];
					String tName = inFields[locusFIdx]+"_"+inFields[targetFIdx];
					int tIdx = getTargetIndex (tName);
					String ref = allTargets[tIdx].getTargetRefSeq();
					int call = SampleCall.getCallFromString(inFields[callFldIdx]);
					String aa = inFields[aaFldIdx];
					String aaNref = inFields[aaNrefFldIdx];
					String nt = inFields[ntFldIdx];
					String ntNref = inFields[ntNrefFldIdx];
					String aaSummary = inFields[countsFldIdx];
					
					SampleCall ntCall = new SampleCall (call, ref, nt, ntNref, null);
					AminoSampleCall aaCall = new AminoSampleCall (call, aa, aaNref, aaSummary, ntCall);
					targetCalls[tIdx][sIdx] = aaCall;
				}
			} finally {
				tif.close();
			}
		}
		return targetCalls;
	}
	
	private SampleTargetResult[][] readAllSampleResults (Sample[] samples) throws AnalysisException {

		SampleTargetResult[][] targetResults = new SampleTargetResult[allTargets.length][samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
				targetResults[tIdx][sIdx] = new SampleTargetResult(allTargets[tIdx], samples[sIdx]);
			}
		}
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
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
		
		int[][] alignedCount    = new int[allTargets.length][samples.length];
		int[][] misalignedCount = new int[allTargets.length][samples.length];
		int[][] coveringCount   = new int[allTargets.length][samples.length];
		int[][] callsCount      = new int[allTargets.length][samples.length];
		int[][] lowQualityCount = new int[allTargets.length][samples.length];
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
			File sampleFile = new File (sampleFolder, sample.getName()+".locusCoverage.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFldIdx      = tif.getFieldIndex("Locus");
			int targetFldIdx     = tif.getFieldIndex("Target");
			int alignedFldIdx    = tif.getFieldIndex("Aligned");
			int misalignedFldIdx = tif.getFieldIndex("Misaligned");
			int coveringFldIdx   = tif.getFieldIndex("Covering");
			int callsFldIdx      = tif.getFieldIndex("Calls");
			int lowQualityFldIdx = tif.getFieldIndex("LowQuality");

			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					//String sampleName = inFields[sampleFIdx];
					String locusName = inFields[locusFldIdx];
					String targetName = inFields[targetFldIdx];
					int tIdx = getTargetIndex (locusName+"_"+targetName);
					alignedCount[tIdx][sIdx]    = Integer.parseInt(inFields[alignedFldIdx]);
					misalignedCount[tIdx][sIdx] = Integer.parseInt(inFields[misalignedFldIdx]);
					coveringCount[tIdx][sIdx]   = Integer.parseInt(inFields[coveringFldIdx]);
					callsCount[tIdx][sIdx]      = Integer.parseInt(inFields[callsFldIdx]);
					lowQualityCount[tIdx][sIdx] = Integer.parseInt(inFields[lowQualityFldIdx]);
				}
			} finally {
				tif.close();
			}
		}
		
		// Write out the results by locus
		for (Locus locus : loci) {
			TableOutput locusOut = new TableOutput (outRootFolder, "LocusCoverage."+locus.getName()+".tab", LOCUS_COUNTS_HEADERS, 64 * 1024);
			TargetLocus tLocus = (TargetLocus)locus;
			Target[] targets = tLocus.getTargets();
			for (Target target : targets) {
				int tIdx = getTargetIndex (locus.getName()+"_"+target.getName());
				for (int sIdx = 0; sIdx < samples.length; sIdx++) {
					locusOut.newRow();
					locusOut.appendValue(samples[sIdx].getBatch());
					locusOut.appendValue(samples[sIdx].getName());
					locusOut.appendValue(locus.getName());
					locusOut.appendValue(target.getName());
					locusOut.appendValue(alignedCount[tIdx][sIdx]);
					locusOut.appendValue(misalignedCount[tIdx][sIdx]);
					locusOut.appendValue(coveringCount[tIdx][sIdx]);
					locusOut.appendValue(callsCount[tIdx][sIdx]);
					locusOut.appendValue(lowQualityCount[tIdx][sIdx]);
				}
			}
			locusOut.close();
		}
	}
	
	private void outputSampleCalls (AminoSampleCall[] targetCalls, Sample[] samples, Target target, File outFolder, String filename) throws AnalysisException {
		int bufferSize = 64 * 1024;
		String[] headers = new String[] {"Batch","Sample","Call","Alleles","NtAlleles","AlleleReads"};
		TableOutput sampleCallsOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < targetCalls.length; sIdx++) {
			AminoSampleCall call = targetCalls[sIdx];
			boolean isLenient = call.isLenient();
			String callString = isLenient ? "["+call.getCallString()+"]" : call.getCallString();
			sampleCallsOut.newRow();
			sampleCallsOut.appendValue(samples[sIdx].getBatch());
			sampleCallsOut.appendValue(samples[sIdx].getName());
			sampleCallsOut.appendValue(callString);
			sampleCallsOut.appendValue(call.getAminoAllele());
			sampleCallsOut.appendValue(call.getAllele());
			sampleCallsOut.appendValue(call.getAminoAlleleSummary());
		}
		sampleCallsOut.close();
	}
	

    private void outputSampleCallsAllTargets (Target[] allTargets, Sample[] samples, AminoSampleCall[][] allCalls) throws AnalysisException {

    	int bufferSize = 1024 * 1024;
    	// Get the file headers
		ArrayList<String> headerList = new ArrayList<String>();
		headerList.add("Batch");
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
			sampleCallsOut.appendValue(samples[sIdx].getBatch());
			sampleCallsOut.appendValue(samples[sIdx].getName());
			
			sampleNrefCallsOut.newRow();
			sampleNrefCallsOut.appendValue(samples[sIdx].getBatch());
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
		String[] headers = new String[aStats.length+2];
		headers[0] = "Batch";
		headers[1] = "Sample";
		for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
			headers[aIdx+2] = aStats[aIdx].displayLabel;
		}
		TableOutput alleleReadsOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			SampleTargetResult sampleResult = sampleResults[sIdx];
			alleleReadsOut.newRow();
			alleleReadsOut.appendValue(sampleResult.getSample().getBatch());
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
		String[] headers = new String[]{"Batch","Sample","Alleles","AlleleCount"};
		TableOutput targetSummaryOut = new TableOutput (outFolder, filename, headers, bufferSize);
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			SampleTargetResult sResult = sampleResults[sIdx];
			targetSummaryOut.newRow();
			targetSummaryOut.appendValue(sResult.getSample().getBatch());
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
			if (args.length < 6) {
				log.error("Usage: org.cggh.bam.codon.CodonAnalysis$SingleSample <configFile> <batchId> <sampleId> <bamFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			String batchId = args[1];					log.info("BatchId: "+batchId);
			String sampleId = args[2];					log.info("SampleId: "+sampleId);
			File sampleBamFile = new File(args[3]);	    log.info("SampleBamFile: "+sampleBamFile.getAbsolutePath());
			File refFastaFile = new File(args[4]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[5]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample sample = new Sample (batchId, sampleId, sampleBamFile);
				CodonAnalysis task = new CodonAnalysis(configFile, refFastaFile, rootFolder);
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
			if (args.length < 4) {
				log.error("Usage: org.cggh.bam.codon.CodonAnalysis$MultiSample <configFile> <sampleListFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[3]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			
			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			
			try {	
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				CodonAnalysis task = new CodonAnalysis(configFile, refFastaFile, rootFolder);
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
			if (args.length < 4) {
				log.error("Usage: org.cggh.bam.codon.CodonAnalysis$MergeResults <configFile> <sampleListFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[3]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				CodonAnalysis task = new CodonAnalysis(configFile, refFastaFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
