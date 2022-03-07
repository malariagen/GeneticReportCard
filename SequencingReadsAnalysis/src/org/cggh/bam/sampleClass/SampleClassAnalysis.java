package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
//import org.cggh.bam.genotyping.*;
import org.cggh.bam.sampleClass.ClassTarget.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class SampleClassAnalysis extends SampleTargetAnalysis  {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private SampleClassConfig     config;
	//private LabelCounterGenotyper genotyper;
	
	private HashMap<String,Integer>[]    alleleIdxTables;
	
	public SampleClassAnalysis (File configFile, File refFastaFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, outRootFolder);

		// Parse configuration file
		config = new SampleClassConfig (configFile);
		log.info(config.getPrintableDisplay());
		registerLoci(config.getLoci());
		
		alleleIdxTables = buildAlleleIndexTables ();
	}

	/* **********************************************************************
	 * Single sample processing
	 * **********************************************************************
	 */
	public void analyzeSample (Sample sample) throws AnalysisException  {
		log.info("Starting " + sample.getName());
		try {
			File outFolder = getSampleSubfolder (outRootFolder, sample, true);
			SampleClassAnalyzer analyzer = new SampleClassAnalyzer (config, outFolder);
			SampleCall call = analyzer.analyzeSample(sample);
			
			// Analyze unlisted alleles: if they are very similar to listed ones, and can be 
			// assigned to a sample class, add them to the relevant counter
			analyzeUnlistedAlleles (call);
			
			// Write out the results
			outputSampleResults (call);
			
		} catch (Exception e) {
			log.info("Aborting " + sample.getName());
			log.error("Error processing BAM file for sample "+ sample.getName() + ": "+e);
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}
	
	private static final String[] LISTED_ALLELES_HEADERS = {"Batch","Sample","Locus","Target","Allele","Count","TargetCall"};
	private static final String[] UNLISTED_ALLELES_HEADERS = {"Batch","Sample","Locus","Target","Allele","Count","Proportion","Closest","Diff"};
	
	/*
	 * Write out the results for this sample into two files: one of counts of listed sample class-specific alleles,
	 * and one for alleles that were not listed.
	 */
	protected void outputSampleResults (SampleCall call) throws AnalysisException, IOException  {
		Sample sample = call.getSample();
		File outFolder = getSampleSubfolder (outRootFolder, sample, true);
		SampleTargetResult[] targetResults = call.getTargetResults();
		
		// Write out the class calls
		String[] callHeaders = TextUtilities.mergeStringLists(new String[]{"Batch","Sample","Class"}, allTargetNames);
		TableOutput callOut = new TableOutput (outFolder, sample.getName()+".classes.tab", callHeaders, 64 * 1024);		
		callOut.newRow();
		callOut.appendValue(sample.getBatch());
		callOut.appendValue(sample.getName());
		callOut.appendValue(call.getClassCall());
		for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
			callOut.appendValue(targetResults[tIdx].getCall());
		}
		callOut.close();
		
		// Write out the allele read counts
		TableOutput alleleSetOut = new TableOutput (outFolder, sample.getName()+".classAlleles.tab", LISTED_ALLELES_HEADERS, 64 * 1024);
		for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
			SampleTargetResult targetResult = targetResults[tIdx];
			ClassTarget target = targetResult.getTarget();
			// Write out allele counts for each allele in the target allele set
			ClassAllele[] alleles = target.getAlleles();
			LabelCounter[] alleleCounters = targetResult.getAlleleCounters();
			String[] classCalls = targetResult.getClassTargetCalls();
			for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
				int alleleReads = alleleCounters[aIdx].getCount();
				alleleSetOut.newRow();
				alleleSetOut.appendValue(sample.getBatch());
				alleleSetOut.appendValue(sample.getName());
				alleleSetOut.appendValue(target.getLocus().getName());
				alleleSetOut.appendValue(target.getName());
				alleleSetOut.appendValue(alleles[aIdx].getName());
				alleleSetOut.appendValue(alleleReads);
				alleleSetOut.appendValue(classCalls[aIdx]);
			}
		}
		alleleSetOut.close();
	}
	

	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	protected void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {

		// Read in all the allele counts and calls for all targets and all samples
		ClassAlleleCounts[] targetCounts = readAllClassAlleleCounts (samples);
		SampleClassCall[] calls = readAllClassCalls (samples);
		
		// Write the calls out
		writeAggregateClassCalls (calls);

		// Write the sample read counts for the class alleles for all loci
		TargetCall[][] targetCalls = new TargetCall[samples.length][allTargets.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			ClassAlleleCounts tCounts = targetCounts[tIdx];
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				targetCalls[sIdx][tIdx] = tCounts.getCall(sIdx);
			}
		}
		writeAggregateSampleCounts (samples, targetCalls);
				
		// Write out the results by target
		int[][] classReadTotals = new int[allTargets.length][];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			ClassAlleleCounts tCounts = targetCounts[tIdx];
			tCounts.writeCounts(outRootFolder);
			classReadTotals[tIdx] = tCounts.getTotalCounts();
		}
		
		// Write out unlisted alleles to file for all samples
		writeUnlistedAllelesByTarget (samples, classReadTotals);
	}
	
	/* **********************************************************************
	 * Reading in and aggregation of target classes
	 * **********************************************************************
	 */
	private SampleClassCall[] readAllClassCalls (Sample[] samples) throws AnalysisException {
		
		SampleClassCall[] results = new SampleClassCall[samples.length];
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
			File sampleFile = new File (sampleFolder, sample.getName()+".classes.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int classFIdx       = tif.getFieldIndex("Class");
			int[] targetFIdxes = new int[allTargets.length]; 
			for (int idx = 0; idx < allTargetNames.length; idx++) {
				targetFIdxes[idx] = tif.getFieldIndex(allTargetNames[idx]);
			}
			try {
				String[] inFields = tif.getNextValidLine();
				if (inFields == null) {
					break;
				}
				String classCall = inFields[classFIdx];
				if ("-".equals(classCall)) {
					classCall = null;
				}
				String[] targetCalls = new String[allTargets.length];
				for (int idx = 0; idx < allTargets.length; idx++) {
				    int fIdx = targetFIdxes[idx];
				    targetCalls[idx] = inFields[fIdx];
					if ("-".equals(targetCalls[idx])) {
						targetCalls[idx] = null;
					}
				}
				results[sIdx] = new SampleClassCall(sample, classCall, targetCalls);
			} finally {
				tif.close();
			}
		}
		return results;
	}

	private class SampleClassCall {
		Sample   sample;
		String   sampleCall;
		String[] targetCalls;

		public SampleClassCall(Sample sample, String sampleCall, String[] targetCalls) {
			super();
			this.sample = sample;
			this.sampleCall = sampleCall;
			this.targetCalls = targetCalls;
		}
	}
	
	private void writeAggregateClassCalls (SampleClassCall[] sampleCalls) throws AnalysisException {
		String[] headers = TextUtilities.mergeStringLists(new String[]{"Batch","Sample","Class"}, allTargetNames);
		TableOutput out = new TableOutput (outRootFolder, "AllSamples-AllTargets.classes.tab", headers, 64 * 1024);
		for (int sIdx = 0; sIdx < sampleCalls.length; sIdx++) {
			out.newRow();
			Sample sample = sampleCalls[sIdx].sample;
			out.appendValue(sample.getBatch());
			out.appendValue(sample.getName());
			out.appendValue(sampleCalls[sIdx].sampleCall);
			for (int tIdx = 0; tIdx < allTargetNames.length; tIdx++) {
				out.appendValue(sampleCalls[sIdx].targetCalls[tIdx]);
			}
		}
		out.close();
	}
	
	/* **********************************************************************
	 * Reading in and aggregation of read counts for target classes
	 * **********************************************************************
	 */
	private ClassAlleleCounts[] readAllClassAlleleCounts (Sample[] samples) throws AnalysisException {
		
		// Create one allele read count table per target, and index the targets by name
		ClassAlleleCounts[] targetCounts = new ClassAlleleCounts[allTargets.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			targetCounts[tIdx] = new ClassAlleleCounts (tIdx, samples);
		}
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
			File sampleFile = new File (sampleFolder, sample.getName()+".classAlleles.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFIdx       = tif.getFieldIndex("Locus");
			int targetFIdx      = tif.getFieldIndex("Target");
			int alleleFldIdx    = tif.getFieldIndex("Allele");
			int countFldIdx     = tif.getFieldIndex("Count");
			int tCallFieldIdx   = tif.getFieldIndex("TargetCall");
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					int count = Integer.parseInt(inFields[countFldIdx]);
					if (count == 0) {
						continue;
					}
					String classCall = inFields[tCallFieldIdx];
					if ("-".equals(classCall)) {
						continue;
					}
					String tName = inFields[locusFIdx]+"_"+inFields[targetFIdx];
					int tIdx = 	getTargetIndex (tName);
					ClassAlleleCounts tCounts = targetCounts[tIdx];
					String allele = inFields[alleleFldIdx];
					tCounts.setCount(sIdx, allele, count);
					tCounts.addClassCall(sIdx, classCall, count);
				}
			} finally {
				tif.close();
			}
		}
		return targetCounts;
	}
	
	private class ClassAlleleCounts {
		Sample[]                samples;
		int                     targetIdx;
		ClassAllele[]           alleles;
		int[][]                 counts;
		TargetCall[]            calls;
		HashMap<String,Integer> alleleIdxTable;
		
		public ClassAlleleCounts (int targetIdx, Sample[] samples) {
			this.samples = samples;
			this.targetIdx = targetIdx;
			this.alleles = ((ClassTarget)allTargets[targetIdx]).getAlleles();
			this.counts = new int[samples.length][alleles.length];
			this.calls = new TargetCall[samples.length];
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				this.calls[sIdx] = new TargetCall();
			}
			this.alleleIdxTable = alleleIdxTables[targetIdx];
		}
		
		public void setCount (int sampleIdx, String allele, int count) {
			int alleleIdx = alleleIdxTable.get(allele);
			counts[sampleIdx][alleleIdx] = count;
		}
		
		public TargetCall getCall(int sampleIdx) {
			return calls[sampleIdx];
		}
		
		public void addClassCall(int sampleIdx, String classCall, int readCount) {
			calls[sampleIdx].addClassCall(classCall, readCount);
		}

		public int[] getTotalCounts () throws AnalysisException {
			int[] totals = new int[samples.length];
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				int total = 0;
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					total += counts[sIdx][aIdx];
				}
				totals[sIdx] = total;
			}
			return totals;
		}
		
		public void writeCounts (File outFolder) throws AnalysisException {
			String[] alleleHeaders = new String[alleles.length+2];
			alleleHeaders[0] = "Sample";
			for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
				ClassAllele allele = alleles[aIdx];
				alleleHeaders[aIdx+1] = allele.getName();
			}
			alleleHeaders[alleles.length+1] = "SampleClass";
			String filename = "AllSamples-"+allTargetNames[targetIdx]+".classAlleles.tab";

			TableOutput alleleSetOut = new TableOutput (outFolder, filename, alleleHeaders, 64 * 1024);
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				alleleSetOut.newRow();
				alleleSetOut.appendValue(samples[sIdx].getName());
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					alleleSetOut.appendValue(counts[sIdx][aIdx]);
				}
				alleleSetOut.appendValue(calls[sIdx].getCallLabel());
			}
			alleleSetOut.close();			
		}
	}
	
	private void writeAggregateSampleCounts (Sample[] samples, TargetCall[][] targetCalls) throws AnalysisException {
		String[] headers = TextUtilities.mergeStringLists(new String[]{"Batch","Sample"}, allTargetNames);
		TableOutput out = new TableOutput (outRootFolder, "AllSamples-AllTargets.counts.tab", headers, 64 * 1024);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			out.newRow();
			out.appendValue(samples[sIdx].getBatch());
			out.appendValue(samples[sIdx].getName());
			for (int tIdx = 0; tIdx < allTargetNames.length; tIdx++) {
				out.appendValue(targetCalls[sIdx][tIdx].getCountsLabel());
			}
		}
		out.close();
	}
	

	private class TargetCall {
		ArrayList<TargetClassCall> callList = new ArrayList<TargetClassCall>();
		
		public void addClassCall (String classCall, int readCount) {
			if (classCall == null) {
				return;
			}
			callList.add(new TargetClassCall(classCall, readCount));
		}
		
		public String getCountsLabel () {
			String callStr = "-";
			if (!callList.isEmpty()) {
				for (int cIdx = 0; cIdx < callList.size(); cIdx++) {
					String cStr = callList.get(cIdx).getCountsLabel();
					callStr = (cIdx == 0) ? cStr : callStr+","+cStr;
				}
			}
			return callStr;
		}
		
		public String getCallLabel () {
			String callStr = "-";
			if (!callList.isEmpty()) {
				for (int cIdx = 0; cIdx < callList.size(); cIdx++) {
					String cStr = callList.get(cIdx).getCallLabel();
					callStr = (cIdx == 0) ? cStr : callStr+","+cStr;
				}
			}
			return callStr;
		}
	}
	
	private class TargetClassCall {
		String targetClass;
		int    readCount;
		
		public TargetClassCall(String targetClass, int readCount) {
			this.targetClass = targetClass;
			this.readCount = readCount;
		}
		public String getCountsLabel () {
			return (Integer.toString(readCount));
		}
		public String getCallLabel () {
			return (targetClass);
		}
	}
	
	private HashMap<String,Integer>[] buildAlleleIndexTables () {
		@SuppressWarnings("unchecked")
		HashMap<String,Integer>[] result = new HashMap[allTargets.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			HashMap<String,Integer> alleleIdxTable = new HashMap<String,Integer>();
			ClassAllele[] alleles = ((ClassTarget)allTargets[tIdx]).getAlleles();
			for (int i = 0; i < alleles.length; i++) {
				alleleIdxTable.put(alleles[i].getName(), i);
			}
			result[tIdx] = alleleIdxTable;
		}
		return result;
	}

	/* ==========================================================
	 * Unlisted alleles analysis
	 * ==========================================================
	 */
	private void analyzeUnlistedAlleles (SampleCall call) throws AnalysisException, IOException  {
		SampleTargetResult[] targetResults = call.getTargetResults();

		UnlistedAllele[][] allUnlistedAlleles = new UnlistedAllele[allTargets.length][];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			SampleTargetResult targetResult = targetResults[tIdx];
			ClassTarget target = targetResult.getTarget();
			LabelCounter[] alleleSetCounters = targetResult.getAlleleCounters();
			double totalReads = targetResult.getTotalReadCount();
			
			// Analyze unlisted alleles: if the allele is very similar to a listed allele (1 mismatch max),
			// and more different from other sample class alleles, and it is in at least 10 reads, then merge it with
			// the class with the most similar allele (add reads to the counter)
			ArrayList<UnlistedAllele> unlistedAlleleList = new ArrayList<UnlistedAllele>();
			LabelCounter[] uCounters = targetResult.getUnlistedAlleleCounters();
			for (int uIdx = 0; uIdx < uCounters.length; uIdx++) {
				String uAllele = uCounters[uIdx].getLabel();
				int uCount = uCounters[uIdx].getCount();
				TargetAlleleSimilarity sim = new TargetAlleleSimilarity(uAllele, target);
				if ((uCount >= 10) && sim.callBySimilarity) {
					alleleSetCounters[sim.mostSimilarAlleleIdx].add(uCount);
				} else {
					double uProp = ((double)uCount) / totalReads;
					unlistedAlleleList.add(new UnlistedAllele (target.getLocus().getName(), target.getName(), uAllele, uCount, uProp, sim));
				}
			}
			allUnlistedAlleles[tIdx] = unlistedAlleleList.toArray(new UnlistedAllele[unlistedAlleleList.size()]);
		}
		outputUnlistedAlleles (call.getSample(), allUnlistedAlleles);
	}

	/*
	 * Write out the results for this sample into two files: one of counts of listed class-specific alleles,
	 * and one for alleles that were not listed.
	 */
	protected void outputUnlistedAlleles (Sample sample, UnlistedAllele[][] allTargetUnlistedAlleles) throws AnalysisException, IOException  {
		File outFolder = getSampleSubfolder (outRootFolder, sample, true);
		TableOutput unlistedAllelesOut = new TableOutput (outFolder, sample.getName()+".unlistedAlleles.tab", UNLISTED_ALLELES_HEADERS, 64 * 1024);
		for (int tIdx = 0; tIdx < allTargetUnlistedAlleles.length; tIdx++) {
			UnlistedAllele[] unlistedAlleles = allTargetUnlistedAlleles[tIdx];
			for (int uIdx = 0; uIdx < unlistedAlleles.length; uIdx++) {
				UnlistedAllele a = unlistedAlleles[uIdx];
				unlistedAllelesOut.newRow();
				unlistedAllelesOut.appendValue(sample.getBatch());
				unlistedAllelesOut.appendValue(sample.getName());
				unlistedAllelesOut.appendValue(a.locusName);
				unlistedAllelesOut.appendValue(a.targetName);
				unlistedAllelesOut.appendValue(a.allele);
				unlistedAllelesOut.appendValue(a.readCount);
				unlistedAllelesOut.appendValue(a.readProportion);
				unlistedAllelesOut.appendValue(a.closestAlleleId);
				unlistedAllelesOut.appendValue(a.closestAlleleMismatches);
			}
	    }
		unlistedAllelesOut.close();
	}
	
	
	private void writeUnlistedAllelesByTarget (Sample[] samples, int[][] classReadTotals) throws AnalysisException {
		
		// Make one output file per target and Index the targets by name
		TableOutput[] unlistedOuts = new TableOutput[allTargets.length];
		HashMap<String,Integer> targetIdxTable = new HashMap<String,Integer>(allTargets.length);
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			unlistedOuts[tIdx] = new TableOutput (outRootFolder, "AllSamples-"+allTargetNames[tIdx]+".unlistedAlleles.tab", UNLISTED_ALLELES_HEADERS, 64 * 1024);
			targetIdxTable.put(allTargetNames[tIdx], tIdx);
		}
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		HashMap<String,UnlistedAlleleStats> statsTable = new HashMap<String,UnlistedAlleleStats>();
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample, true);
			File sampleFile = new File (sampleFolder, sample.getName()+".unlistedAlleles.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFIdx   = tif.getFieldIndex("Locus");
			int targetFIdx  = tif.getFieldIndex("Target");
			int alleleFIdx  = tif.getFieldIndex("Allele");
			int countFIdx   = tif.getFieldIndex("Count");
			int propFIdx    = tif.getFieldIndex("Proportion");
			int closestFIdx = tif.getFieldIndex("Closest");
			int diffFIdx    = tif.getFieldIndex("Diff");
			
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					
					// Copy the data into an al-sample aggregated file
					String locusName = inFields[locusFIdx];
					String targetName = inFields[targetFIdx];
					String tName = locusName+"_"+targetName;
					int tIdx = targetIdxTable.get(tName);
					
					unlistedOuts[tIdx].newRow();
					String[] currFields = tif.getCurrentLineFields();
					for (int fIdx = 1; fIdx < currFields.length; fIdx++) {
						unlistedOuts[tIdx].appendValue(currFields[fIdx]);
					}
					
					// Update the statistics for this target allele
					String allele = inFields[alleleFIdx];
					String alleleId = tName+"/"+allele;
					int readCount = Integer.parseInt(inFields[countFIdx]);
					double prop = Double.parseDouble(inFields[propFIdx]);
					
					UnlistedAlleleStats stats = statsTable.get(alleleId);
					if (stats == null) {
						String closestAlleleId = inFields[closestFIdx];
						int diff = Integer.parseInt(inFields[diffFIdx]);
						stats = new UnlistedAlleleStats(locusName, targetName, allele, closestAlleleId, diff);
						statsTable.put(alleleId, stats);
					}
					stats.sampleCount++;
					stats.maxReadCount = (readCount > stats.maxReadCount) ? readCount : stats.maxReadCount;
					stats.maxReadProportion = (prop > stats.maxReadProportion) ? prop : stats.maxReadProportion;
				}
			} finally {
				tif.close();
				for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
					unlistedOuts[tIdx].close();
				}
			}
		}
		
		UnlistedAlleleStats[] allStats = statsTable.values().toArray(new UnlistedAlleleStats[statsTable.values().size()]);
		Arrays.sort(allStats);
		String[] headers = {"Locus","Target","Allele","SampleCount","MaxReadProp","MaxReadCount","Closest","Diff"};
		String filename = "AllSamples-UnlistedAllelesStats.tab";
		TableOutput statOut = new TableOutput (outRootFolder, filename, headers, 64 * 1024);
		for (int sIdx = 0; sIdx < allStats.length; sIdx++) {
			if (allStats[sIdx].maxReadCount == 1) {
				continue;  // Trim out singletons, at least
			}
			statOut.newRow();
			statOut.appendValue(allStats[sIdx].locusName);
			statOut.appendValue(allStats[sIdx].targetName);
			statOut.appendValue(allStats[sIdx].allele);
			statOut.appendValue(allStats[sIdx].sampleCount);
			statOut.appendValue(allStats[sIdx].maxReadProportion);
			statOut.appendValue(allStats[sIdx].maxReadCount);
			statOut.appendValue(allStats[sIdx].closestAlleleId);
			statOut.appendValue(allStats[sIdx].closestAlleleMismatches);
		}
		statOut.close();			
	}

	
	private class UnlistedAllele implements Comparable<UnlistedAllele>{
		String locusName;
		String targetName;
		String allele;
		String closestAlleleId;
		int    closestAlleleMismatches;
		int    readCount;
		double readProportion;


		public UnlistedAllele (String locusName, String targetName, String allele, int readCount, double readProportion, TargetAlleleSimilarity sim) {
			this.locusName = locusName;
			this.targetName = targetName;
			this.allele = allele;
			this.readCount = readCount;
			this.readProportion = readProportion;
			this.closestAlleleId = sim.getMostSimilarAllele().getName();
			this.closestAlleleMismatches = sim.leastMismatches;
		}

		@Override
		public int compareTo(UnlistedAllele o) {
			int result = locusName.compareTo(o.locusName);
			if (result != 0) return result;
			result = targetName.compareTo(o.targetName);
			if (result != 0) return result;
			return( o.readCount - readCount);
		}
	}
	
	
	private class TargetAlleleSimilarity {
		ClassTarget target;
		int mostSimilarAlleleIdx;
		int leastMismatches;
		boolean callBySimilarity;
		
		public TargetAlleleSimilarity(String allele, ClassTarget target) {
			this.target = target;
			char[] seq = allele.toCharArray();
			int len = seq.length;
			
			mostSimilarAlleleIdx = -1;
			leastMismatches = 1000;
			int secondLeastMismatches = 1000;
			
			ClassAllele[] alleles = target.getAlleles();
			for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
				int alleleMismatches = 1000;
				String[] aSeqs = alleles[aIdx].getSequences();
				for (int sIdx = 0; sIdx < aSeqs.length; sIdx++) {
					char[] testSeq = aSeqs[sIdx].toCharArray();
					int seqMismatches = 0;
					for (int i = 0; i < len; i++) {
						if (testSeq[i] != seq[i]) {
							seqMismatches++;
						}
					}
					if (seqMismatches < alleleMismatches) {
						alleleMismatches = seqMismatches;
					}
				}
				if (alleleMismatches < leastMismatches) {
					secondLeastMismatches = leastMismatches;
					leastMismatches = alleleMismatches;
					mostSimilarAlleleIdx = aIdx;
				} else if (alleleMismatches < secondLeastMismatches) {
					secondLeastMismatches = alleleMismatches;
				}
			}
			
			// If there is a single mismatch from the most similar allele, and the second most similar has at least 3, call the allele
			callBySimilarity = ((leastMismatches == 1) && (secondLeastMismatches >= 3));
		}
		
		public ClassAllele getMostSimilarAllele () {
			return target.getAlleles()[mostSimilarAlleleIdx];
		}
	}
	
	private class UnlistedAlleleStats implements Comparable<UnlistedAlleleStats>{
		String locusName;
		String targetName;
		String allele;
		String closestAlleleId;
		int    closestAlleleMismatches;
		int    sampleCount;
		int    maxReadCount;
		double maxReadProportion;
		
		public UnlistedAlleleStats (String locusName, String targetName, String allele, String closestAlleleId, int closestAlleleMismatches) {
			this.locusName = locusName;
			this.targetName = targetName;
			this.allele = allele;
			this.closestAlleleId = closestAlleleId;
			this.closestAlleleMismatches = closestAlleleMismatches;
		}
				
		@Override
		public int compareTo(UnlistedAlleleStats o) {
			int result = locusName.compareTo(o.locusName);
			if (result != 0) return result;
			result = targetName.compareTo(o.targetName);
			if (result != 0) return result;
			result = o.maxReadCount - maxReadCount;
			if (result != 0) return result;
			return (o.sampleCount - sampleCount);
		}
	}
	


	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		
		public static void main(String[] args) {
			if (args.length < 6) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$SingleSample <configFile> <batchId> <sampleId> <bamFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			String batchId = args[1];					log.info("BatchId: "+batchId);
			String sampleId = args[2];					log.info("SampleId: "+sampleId);
			File sampleBamFile = new File(args[3]);	    log.info("SampleBamFile: "+sampleBamFile.getAbsolutePath());
			File refFastaFile = new File(args[4]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[5]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				SampleClassAnalysis task = new SampleClassAnalysis(configFile, refFastaFile, rootFolder);
				Sample sample = new Sample (batchId, sampleId, sampleBamFile);
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
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MultiSample <configFile> <sampleListFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[3]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			//maxThreads = 1;
			
			try {	
				SampleClassAnalysis task = new SampleClassAnalysis(configFile, refFastaFile, rootFolder);
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				multi.execute(task);	// Calls analyzeSample() for each sample
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
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MergeResults <configFile> <sampleListFile> <refFasta> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File rootFolder = new File(args[3]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				SampleClassAnalysis task = new SampleClassAnalysis(configFile, refFastaFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
