package org.cggh.bam.species;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.bam.target.alleleClasses.*;
import org.cggh.bam.target.alleleClasses.SampleAlleleClassAnalyzer.*;
import org.cggh.bam.target.alleleClasses.AlleleClassTarget.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class SpeciesAnalysis extends SampleTargetAnalysis  {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private SpeciesConfig            config;
	
	public SpeciesAnalysis (File configFile, File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, chrMapFile, outRootFolder);

		// Parse configuration file
		config = new SpeciesConfig (configFile);
		registerLoci(config.getLoci());
	}

	/* **********************************************************************
	 * Single sample processing
	 * **********************************************************************
	 */
	public void analyzeSample (Sample sample) throws AnalysisException  {
		log.info("Starting " + sample.getName());
		try {
			SampleAlleleClassAnalyzer analyzer = new SampleAlleleClassAnalyzer ((AlleleClassLocus[])loci, sample);
			SampleAlleleClassResults sr = analyzer.analyzeSample();
			
			// Analyze unlisted alleles: if they are very similar to listed ones, and can be 
			// assigned to a species, add them to the relevant counter
			analyzeUnlistedAlleles (sr);
			
			// Write out the results
			outputSampleResults (sr);
			
		} catch (Exception e) {
			log.info("Aborting " + sample.getName());
			log.error("Error processing BAM file for sample "+ sample.getName() + ": "+e);
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}
	
	private static final String[] LISTED_ALLELES_HEADERS = {"Sample","Locus","Target","Allele","Count","Species"};
	private static final String[] UNLISTED_ALLELES_HEADERS = {"Sample","Locus","Target","Allele","Count","Proportion","Closest","Diff"};
	
	/*
	 * Write out the results for this sample into two files: one of counts of listed species-specific alleles,
	 * and one for alleles that were not listed.
	 */
	protected void outputSampleResults (SampleAlleleClassResults sr) throws AnalysisException, IOException  {
		Sample sample = sr.getSample();
		File outFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);

		TableOutput alleleSetOut = new TableOutput (outFolder, sample.getName()+".speciesAlleles.tab", LISTED_ALLELES_HEADERS, 64 * 1024);		

		SampleAlleleClassLocusResult[] locusResults = sr.getLocusResults();
		for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
			SampleAlleleClassLocusResult locusResult = locusResults[lIdx];
			AlleleClassLocus locus = locusResult.getLocus();
			
			SampleAlleleClassResult[] targetResults = locusResult.getTargetResults();
			for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
				SampleAlleleClassResult targetResult = targetResults[tIdx];
				AlleleClassTarget target = targetResult.getTarget();
				
				// Write out allele counts for each allele in the target allele set
				ClassAllele[] alleles = target.getAlleles();
				LabelCounter[] alleleCounters = targetResult.getAlleleSetCounters();
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					alleleSetOut.newRow();
					alleleSetOut.appendValue(sample.getName());
					alleleSetOut.appendValue(locus.getName());
					alleleSetOut.appendValue(target.getName());
					
					String alleleSeq = alleles[aIdx].getName();
					int count = alleleCounters[aIdx].getCount();
					String speciesStr = (count > 1) ? alleleCounters[aIdx].getLabel() : "-";
					alleleSetOut.appendValue(alleleSeq);
					alleleSetOut.appendValue(count);
					alleleSetOut.appendValue(speciesStr);
				}
			}
		}
		alleleSetOut.close();
	}
	

	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	protected void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {

		// Go through all the samples, reading in all the allele counts for each target
		ClassAlleleCounts[] targetCounts = readAllClassAlleleCounts (samples);
		
		// Write out the results by target
		int[][] classReadTotals = new int[allTargets.length][];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			ClassAlleleCounts tCounts = targetCounts[tIdx];
			tCounts.writeCounts(outRootFolder);
			classReadTotals[tIdx] = tCounts.getTotalCounts();
		}
		
		// Write out unlisted alleles to file for all samples
		writeUnlistedAllelesByTarget (samples, classReadTotals);
		
		// Get the species alleles for all loci
		String[][] speciesAlleles = new String[samples.length][allTargets.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			ClassAlleleCounts tCounts = targetCounts[tIdx];
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				speciesAlleles[sIdx][tIdx] = tCounts.species[sIdx];
			}
		}

		// Process one class at a time to get the overall call for the sample
		String[] sampleSpecies = callSamples (samples, (AlleleClassTarget[])allTargets, speciesAlleles);

		// Write the calls out
		String[] speciesHeaders = TextUtilities.mergeStringLists(new String[] {"Sample", "Species", "Consensus"}, allTargetNames);
		TableOutput speciesOut = new TableOutput (outRootFolder, "AllSamples-AllTargets.species.tab", speciesHeaders, 64 * 1024);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			speciesOut.newRow();
			speciesOut.appendValue(samples[sIdx].getName());
			String call = sampleSpecies[sIdx];
			String cons;
			if (call != null) {
				cons = call.contains("*") ? "N" : "Y";
			} else {
				call = "-";
				cons = "-";
			}
			speciesOut.appendValue(call);
			speciesOut.appendValue(cons);
			for (int tIdx = 0; tIdx < allTargetNames.length; tIdx++) {
				speciesOut.appendValue(speciesAlleles[sIdx][tIdx]);
			}
		}
		speciesOut.close();
	}
	
	
	private ClassAlleleCounts[] readAllClassAlleleCounts (Sample[] samples) throws AnalysisException {
		
		// Create one allele read count table per target, and index the targets by name
		ClassAlleleCounts[] targetCounts = new ClassAlleleCounts[allTargets.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			targetCounts[tIdx] = new ClassAlleleCounts (tIdx, samples);
		}
		
		// Go through all the samples, reading in all the allele counts for each target, and put them in the allele read count objects
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
			File sampleFile = new File (sampleFolder, sample.getName()+".speciesAlleles.tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			TableInput tif = new TableInput (sampleFile);
			int locusFIdx       = tif.getFieldIndex("Locus");
			int targetFIdx      = tif.getFieldIndex("Target");
			int alleleFldIdx    = tif.getFieldIndex("Allele");
			int countFldIdx     = tif.getFieldIndex("Count");
			int speciesFieldIdx = tif.getFieldIndex("Species");
			try {
				while (true) {
					String[] inFields = tif.getNextValidLine();
					if (inFields == null) {
						break;
					}
					//String sampleName = inFields[sampleFIdx];
					String tName = inFields[locusFIdx]+"_"+inFields[targetFIdx];
					int tIdx = 	getTargetIndex (tName);
					String allele = inFields[alleleFldIdx];
					int count = Integer.parseInt(inFields[countFldIdx]);
					String speciesCall = inFields[speciesFieldIdx];
					
					ClassAlleleCounts tCounts = targetCounts[tIdx];
					tCounts.setCount(sIdx, allele, count);
					tCounts.addSpeciesCall(sIdx, "-".equals(speciesCall) ? null : speciesCall);
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
		String[]                species;
		HashMap<String,Integer> alleleIdxTable;
		
		public ClassAlleleCounts (int targetIdx, Sample[] samples) {
			this.samples = samples;
			this.targetIdx = targetIdx;
			this.alleles = ((AlleleClassTarget)allTargets[targetIdx]).getAlleles();
			this.counts = new int[samples.length][alleles.length];
			this.species = new String[samples.length];
			alleleIdxTable = new HashMap<String,Integer>(alleles.length);
			for (int i = 0; i < alleles.length; i++) {
				alleleIdxTable.put(alleles[i].getName(), i);
			}
		}
		
		public void setCount (int sampleIdx, String allele, int count) {
			int alleleIdx = alleleIdxTable.get(allele);
			counts[sampleIdx][alleleIdx] = count;
		}
		
		public void addSpeciesCall (int sampleIdx, String speciesCall) {
			if (speciesCall == null) {
				return;
			}
			String speciesStr = species[sampleIdx];
			if (speciesStr == null) {
				speciesStr = speciesCall;
			} else {
				speciesStr += ","+speciesCall;
			}
			species[sampleIdx] = speciesStr;
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
			alleleHeaders[alleles.length+1] = "Species";
			String filename = "AllSamples-"+allTargetNames[targetIdx]+".speciesAlleles.tab";

			TableOutput alleleSetOut = new TableOutput (outFolder, filename, alleleHeaders, 64 * 1024);
			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				alleleSetOut.newRow();
				alleleSetOut.appendValue(samples[sIdx].getName());
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					alleleSetOut.appendValue(counts[sIdx][aIdx]);
				}
				alleleSetOut.appendValue(species[sIdx]);
			}
			alleleSetOut.close();			
		}
	}
	
	/* ==========================================================
	 * Unlisted alleles analysis
	 * ==========================================================
	 */
	private void analyzeUnlistedAlleles (SampleAlleleClassResults sr) throws AnalysisException, IOException  {

		UnlistedAllele[][] allUnlistedAlleles = new UnlistedAllele[allTargets.length][];
		int allTargetIdx = 0;
		SampleAlleleClassLocusResult[] locusResults = sr.getLocusResults();
		for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
			SampleAlleleClassLocusResult locusResult = locusResults[lIdx];
			AlleleClassLocus locus = locusResult.getLocus();
			SampleAlleleClassResult[] targetResults = locusResult.getTargetResults();
			for (int tIdx = 0; tIdx < targetResults.length; tIdx++) {
				SampleAlleleClassResult targetResult = targetResults[tIdx];
				AlleleClassTarget target = targetResult.getTarget();
				LabelCounter[] alleleSetCounters = targetResult.getAlleleSetCounters();
				double totalReads = targetResult.getTotalReadCount();
				
				// Analyze unlisted alleles: if the allele is very similar to a listed allele (1 mismatch max),
				// and more different from other species alleles, and it is in at least 10 reads, then merge it with
				// the species with the most similar allele (add reads to the counter)
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
						unlistedAlleleList.add(new UnlistedAllele (locus.getName(), target.getName(), uAllele, uCount, uProp, sim));
					}
				}
				allUnlistedAlleles[allTargetIdx++] = unlistedAlleleList.toArray(new UnlistedAllele[unlistedAlleleList.size()]);
			}
		}
		outputUnlistedAlleles (sr.getSample(), allUnlistedAlleles);
	}

	/*
	 * Write out the results for this sample into two files: one of counts of listed species-specific alleles,
	 * and one for alleles that were not listed.
	 */
	protected void outputUnlistedAlleles (Sample sample, UnlistedAllele[][] allTargetUnlistedAlleles) throws AnalysisException, IOException  {
		File outFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
		TableOutput unlistedAllelesOut = new TableOutput (outFolder, sample.getName()+".unlistedAlleles.tab", UNLISTED_ALLELES_HEADERS, 64 * 1024);
		for (int tIdx = 0; tIdx < allTargetUnlistedAlleles.length; tIdx++) {
			UnlistedAllele[] unlistedAlleles = allTargetUnlistedAlleles[tIdx];
			for (int uIdx = 0; uIdx < unlistedAlleles.length; uIdx++) {
				UnlistedAllele a = unlistedAlleles[uIdx];
				unlistedAllelesOut.newRow();
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
			File sampleFolder = getSampleSubfolder (outRootFolder, sample.getName(), true);
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
		AlleleClassTarget target;
		int mostSimilarAlleleIdx;
		int leastMismatches;
		boolean callBySimilarity;
		
		public TargetAlleleSimilarity(String allele, AlleleClassTarget target) {
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
	 * Species Calling from all targets
	 * ==========================================================
	 */
	private String[] callSamples (Sample[] samples, AlleleClassTarget[] allTargets, String[][] speciesAlleles) {
		
		// Get the list of species we're trying to find
		String[] classes = config.getClasses();
		HashMap<String,Integer> classesIdxTable = new HashMap<String,Integer>();
		for (int cIdx = 0; cIdx < classes.length; cIdx++) {
			classesIdxTable.put(classes[cIdx], cIdx);
		}
		
		// Process one class at a time to get the overall call for the sample
		String[] sampleSpecies = new String[samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			int[] specificAlleleCounts = new int[classes.length];
			int[] promiscuousAllelesCounts = new int[classes.length];
			int validCount = 0;
			for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
				String speciesValue = speciesAlleles[sIdx][tIdx];
				if ((speciesValue == null) || speciesValue.equals("-")) {
					continue;
				}
				validCount++;
				String[] alleleLabels = speciesValue.split(",");
				for (String alleleLabel : alleleLabels) {
					String[] alleleClasses = alleleLabel.split("\\|");
					boolean isPromiscuous = alleleClasses.length > 1;
					for (String alleleClass : alleleClasses) {
						Integer idxObj = classesIdxTable.get(alleleClass);
						if (idxObj == null) {
							log.error("Allele '"+alleleClass+"' found at "+allTargets[tIdx]+" in sample "+ samples[sIdx].getName()+ " is not a valid class");
						} else {
							int idx = idxObj.intValue();
							if (isPromiscuous) {
								promiscuousAllelesCounts[idx]++;
							} else {
								specificAlleleCounts[idx]++;
							}
						}
					}
				}
			}
			
			// Do the calling
			String call = null;
			for (int cIdx = 0; cIdx < classes.length; cIdx++) {
				boolean callClass = specificAlleleCounts[cIdx] > 0;
				if (callClass) {
					call = (call == null) ? classes[cIdx] : call+","+classes[cIdx];
					boolean consensus = (specificAlleleCounts[cIdx] + promiscuousAllelesCounts[cIdx]) == validCount;
					if (!consensus) {
						call += "*";
					}
				}
			}
			sampleSpecies[sIdx] = (call == null) ? "-" : call;
		}
		return sampleSpecies;
	}
	
	

	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		
		public static void main(String[] args) {
			if (args.length < 6) {
				log.error("Usage: org.cggh.bam.species.SpeciesAnalysis$SingleSample <configFile> <sampleName> <bamFile> <chrMap> <refFasta> <chrMapFile> <rootFolder>");
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
				SpeciesAnalysis task = new SpeciesAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
				Sample sample = new Sample (sampleId, sampleBamFile, sampleChrMapName);
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
				log.error("Usage: org.cggh.bam.species.SpeciesAnalysis$MultiSample <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);		log.info("ChrMapFile: "+chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			
			try {	
				SpeciesAnalysis task = new SpeciesAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
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
			if (args.length < 5) {
				log.error("Usage: org.cggh.bam.species.SpeciesAnalysis$MergeResults <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		log.info("ConfigFile: "+configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);		log.info("RefFastaFile: "+refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);		log.info("ChrMapFile: "+chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				SpeciesAnalysis task = new SpeciesAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
