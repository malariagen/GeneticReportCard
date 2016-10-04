package org.cggh.bam.genotyping;

import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class TargetsAlignmentAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	private int bufferSize =  4 * 1024 * 1024;
	
	private Sample[]      samples;
	private Locus[] locusConfigs;
	
	private File rootFolder;
	private File lociFolder;
	private File targetsFolder;

	
	
	/* ==========================================================
	 * Samples and SNP lists initialization
	 * ==========================================================
	 */
	public TargetsAlignmentAnalysis (File configFile, File sampleListFile, File refFastaFile, File rootFolder) throws AnalysisException  {

		// Load up the reference genome sequences
		ReferenceGenome.initialize(refFastaFile);
		
		GenotypingConfig config = new GenotypingConfig (configFile);
		locusConfigs = config.getLoci();
		
		// Get the list of samples
		samples = readSamples(sampleListFile);

		// Create the root output directory
		this.rootFolder = rootFolder;
		this.lociFolder    = FileUtilities.checkFolder(rootFolder, "loci", false);
		this.targetsFolder = FileUtilities.checkFolder(rootFolder, "targets", true);
	}


	public void execute () throws AnalysisException, IOException  {
		
		//SampleLocusResult[][] slResults = new SampleLocusResult[samples.length][locusConfigs.length];
		LocusResult[] locusResults = new LocusResult[locusConfigs.length];

		for (int lIdx = 0; lIdx < locusConfigs.length; lIdx++) {
			LocusResult locusResult = new LocusResult(locusConfigs[lIdx], samples);
			log.info("Analyzing locus "+locusConfigs[lIdx].getName());
			ReadsAlignmentRetriever raRetriever = new ReadsAlignmentRetriever (locusConfigs[lIdx], lociFolder);
			Target[] targets = locusConfigs[lIdx].getTargets();
			
			AlignmentTarget[] emptyTargets = new AlignmentTarget[targets.length];
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				emptyTargets[tIdx] = new AlignmentTarget(targets[tIdx]);
			}

			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				
				ReadsAlignment ra = null;
				try {
					ra = raRetriever.retrieveReadsAlignment(samples[sIdx]);					
				} catch (Exception e) {
					log.error ("Error reading reads alignment file for sample "+samples[sIdx].getName()+": "+e);
					ra = null;
				}
				if (ra == null) {
					locusResult.sampleResults[sIdx] = new SampleLocusResult(samples[sIdx], 0, 0, emptyTargets);
					continue;
				}
				
				MappedRead[] mappedReads = ra.getMappedReads();
				int readCount = ra.getReadCount();
				
				AlignmentTarget[] alignTargets = new AlignmentTarget[targets.length];
				for (int tIdx = 0; tIdx < targets.length; tIdx++) {
					Target t = targets[tIdx];
					GenomeRegion tRegion = t.getTargetRegion();
					int tStartPos = tRegion.getStartPos();
					int tEndPos = tRegion.getStopPos();
					int tLen = 1 + tEndPos - tStartPos;

					AlignmentTarget at = new AlignmentTarget(t);
					for (int rIdx = 0; rIdx < readCount; rIdx++) {
						MappedRead r = mappedReads[rIdx];
						int rStartPos = r.getStartPos();
						int rEndPos = rStartPos + r.getSequence().length() - 1;
						if ((rStartPos <= tStartPos) && (rEndPos >= tEndPos)) {
							int tStartOffset = tStartPos - rStartPos;
							String tSeq = r.getSequence().substring(tStartOffset, tStartOffset+tLen);
							if (!tSeq.contains("N")) {
								int minQ = 1000;
								for (int j = 0; j < tLen; j++) {
									int q = r.getPhredScore(tStartOffset+j);
									if (q < minQ) {
										minQ = q;
									}
								}
								
								// If the gene is negative-strand, reverse the sequence
								if (t.isReverse()) {
									tSeq = SequenceUtilities.getReverseComplementSequence(tSeq);					
								}

								// Catalogue the target alleles
								at.alleleCount++;
								at.rawAlleleCounters.increment(tSeq);
								if (minQ >= MappedRead.MIN_PHRED_SCORE) {
									at.alleleCounters.increment(tSeq);
								}
							}
						}
					}
					at.cleanupTargetAlleles();
					alignTargets[tIdx] = at;
				}
				locusResult.sampleResults[sIdx] = new SampleLocusResult(samples[sIdx], readCount, ra.getMisalignedCount(), alignTargets);
			}
			
			// Count the samples per allele in the master lists
			for (int tIdx = 0; tIdx < locusResult.getTargetCount(); tIdx++) {
				TargetAlleleCounters tac = locusResult.alleleSampleCounters[tIdx];
				for (int sIdx = 0; sIdx < samples.length; sIdx++) {
					SampleLocusResult slr = locusResult.sampleResults[sIdx];
					LabelCounter[] ac = slr.alignmentTargets[tIdx].alleleCounters.getSortedCounters();
					for (int aIdx = 0; aIdx < ac.length; aIdx++) {
						tac.increment(ac[aIdx].getLabel());						
					}
				}
			}
			
			// Output the sample summary for this locus
			outputSampleSummaries (locusResult);

		    // Now analyze the target alleles
			SampleLocusResult[] sResults = locusResult.sampleResults;
			for (int tIdx = 0; tIdx < locusResult.getTargetCount(); tIdx++) {
				
				// Compute some target Allele statistics
				AlleleStats[] aStats = computeTargetStats (locusResult, tIdx);
				
				// Now go through and filter the alleles to reject likely artefact
				// The test is MinAlt: homozygous or >=10 reads in at least one sample
				ArrayList<AlleleStats> filteredAllelelStatList = new ArrayList<AlleleStats>();
				for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
					if ((aStats[aIdx].maxReadsFraction == 1.0) || (aStats[aIdx].maxReads >= 10)) {
						filteredAllelelStatList.add(aStats[aIdx]);
					}
				}
				AlleleStats[] filteredAlleleStats = filteredAllelelStatList.toArray(new AlleleStats[filteredAllelelStatList.size()]);
				String[] filteredAlleles = new String[filteredAlleleStats.length];
				for (int aIdx = 0; aIdx < filteredAlleleStats.length; aIdx++) {
					filteredAlleles[aIdx] = filteredAlleleStats[aIdx].allele;
				}
				
				// Go through the samples, keeping only the alleles that passed the MinAlt filter.
				LabelCounterFilter alleleListFilter = new LabelCounterFilterByLabelList (filteredAlleles);
				for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
					sResults[sIdx].alignmentTargets[tIdx].filterAlleleCounters(alleleListFilter);
				}

				// Final step: make a call 
				for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
					sResults[sIdx].alignmentTargets[tIdx].callGenotypes();
				}
				
				// Write out a final allele read count table
				outputAlleleReadCounts (locusResult, tIdx, filteredAlleleStats, "-final");
				
				// Write out an allele stats table
				outputAlleleStats (locusResult, tIdx, filteredAlleleStats, "-final");
				
				// Write out the final calls
				outputSampleCalls (locusResult, tIdx);
			}
			locusResults[lIdx] = locusResult;
		}
		
		// Write out the calls for all targets in one file
		outputSampleCallsAllTargets (locusResults);
	}

	
	
	/* *************************************************************************
	 * Initialization routines
	 * *************************************************************************
	 */
	private Sample[] readSamples(File sampleListFile) throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));
		String[] colNames = new String[]{"Sample"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample(values[0], null);
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}

	/* *************************************************************************
	 * 
	 * *************************************************************************
	 */
	private AlleleStats[] computeTargetStats (LocusResult locusResult, int targetIdx) throws AnalysisException {
		TargetAlleleCounters alleleSampleCounters = locusResult.alleleSampleCounters[targetIdx];
		SampleLocusResult[] sampleResults = locusResult.sampleResults;

		LabelCounter[] aSampleCounters = alleleSampleCounters.getSortedCounters();
		AlleleStats[] aStats = new AlleleStats[aSampleCounters.length];
		for (int aIdx = 0; aIdx < aSampleCounters.length; aIdx++) {
			TargetAlleleCounter ac = (TargetAlleleCounter)aSampleCounters[aIdx];
			aStats[aIdx] = new AlleleStats(ac.getLabel(), ac.fullLabel);
		}
		
		for (int sIdx = 0; sIdx < sampleResults.length; sIdx++) {
			LabelCounters aCounters = sampleResults[sIdx].alignmentTargets[targetIdx].alleleCounters;
			double totalSampleReads = aCounters.getTotal();
			for (int aIdx = 0; aIdx < aSampleCounters.length; aIdx++) {
				AlleleStats stat = aStats[aIdx];
				LabelCounter c = aCounters.getCounter(stat.allele);
				if (c != null) {
					int readCount = c.getCount();
					double readsFraction = ((double)readCount) / totalSampleReads;
					stat.sampleCount++;
					stat.maxReads = (readCount > stat.maxReads) ? readCount : stat.maxReads;
					stat.maxReadsFraction = (readsFraction > stat.maxReadsFraction) ? readsFraction : stat.maxReadsFraction;
				}
			}
		}
		
		// Write out an allele read count table
		outputAlleleReadCounts (locusResult, targetIdx, aStats, "-beforeMinAlt");

		// Write out an allele stats table
		outputAlleleStats (locusResult, targetIdx, aStats, "-beforeMinAlt");

		return aStats;
	}
		
	public class AlleleStats {
		String allele;
		String displayLabel;
		int    sampleCount;
		int    maxReads;
		double maxReadsFraction;
		
		public AlleleStats(String allele, String displayLabel) {
			super();
			this.allele = allele;
			this.displayLabel = displayLabel;
		}
	}

	/* ==========================================================
	 * 
	 * ==========================================================
	 */
	private void outputSampleSummaries (LocusResult locusResult) throws AnalysisException {

		SampleLocusResult[] sResults = locusResult.sampleResults;

		// Write the locus sample summary file
		File outLocusFolder = FileUtilities.checkFolder(lociFolder, locusResult.locus.getName(), true);
		String[] colNames = {"Sample","TotalReads","MisalignedReads"};
		TableOutput locusSummaryOut = new TableOutput (outLocusFolder, "LocusSummaryBySample-beforeMinAlt.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
			SampleLocusResult sResult = sResults[sIdx];
			locusSummaryOut.newRow();
			locusSummaryOut.appendValue(sResult.sample.getName());
			locusSummaryOut.appendValue(sResult.readCount);
			locusSummaryOut.appendValue(sResult.misalignCount);
		}
		locusSummaryOut.close();
		
		// Write the target sample summary file
		for (int tIdx = 0; tIdx < locusResult.getTargetCount(); tIdx++) {
			Target t = locusResult.locus.getTargets()[tIdx];
			File outTargetFolder = FileUtilities.checkFolder(targetsFolder, t.getName(), true);
			colNames = new String[]{"Sample","Alleles","UnfilteredAlleles","AlleleCount"};
			TableOutput targetSummaryOut = new TableOutput (outTargetFolder, "TargetSummaryBySample-beforeMinAlt.tab", colNames, bufferSize);
			for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
				SampleLocusResult sResult = sResults[sIdx];
				targetSummaryOut.newRow();
				targetSummaryOut.appendValue(sResult.sample.getName());
				targetSummaryOut.appendValue(sResult.alignmentTargets[tIdx].alleleCounters.getSummary());
				targetSummaryOut.appendValue(sResult.alignmentTargets[tIdx].rawAlleleCounters.getSummary());
				targetSummaryOut.appendValue(sResult.alignmentTargets[tIdx].alleleCount);
			}
			targetSummaryOut.close();
		}
	}
	
	private void outputAlleleReadCounts (LocusResult locusResult, int targetIdx, AlleleStats[] aStats, String suffix) throws AnalysisException {

		SampleLocusResult[] sResults = locusResult.sampleResults;

		// Write out an allele read count table
		String[] aDisplayLabels = new String[aStats.length];
		for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
			aDisplayLabels[aIdx] = aStats[aIdx].displayLabel;
		}
		
		Target t = locusResult.locus.getTargets()[targetIdx];
		File outTargetFolder = FileUtilities.checkFolder(targetsFolder, t.getName(), true);
		TableOutput alleleReadsOut = new TableOutput (outTargetFolder, "AlleleSampleCount"+suffix+".tab",
				               TextUtilities.mergeStringLists(new String[] {"Sample"}, aDisplayLabels), bufferSize);
		for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
			SampleLocusResult sResult = sResults[sIdx];
			alleleReadsOut.newRow();
			alleleReadsOut.appendValue(sResult.sample.getName());
			for (int aIdx = 0; aIdx < aStats.length; aIdx++) {
				LabelCounter c = sResult.alignmentTargets[targetIdx].alleleCounters.getCounter(aStats[aIdx].allele);
				alleleReadsOut.appendValue((c != null) ? c.getCount() : 0);
			}
		}
		alleleReadsOut.close();
	}

	private void outputAlleleStats (LocusResult locusResult, int targetIdx, AlleleStats[] hStats, String suffix) throws AnalysisException {
		Target t = locusResult.locus.getTargets()[targetIdx];
		File outTargetFolder = FileUtilities.checkFolder(targetsFolder, t.getName(), true);
		TableOutput alleleStatsOut = new TableOutput (outTargetFolder, "AlleleStats"+suffix+".tab", 
                new String[] {"Allele","SampleCount","MaxReads","MaxReadFraction"}, bufferSize);
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			alleleStatsOut.newRow();
			alleleStatsOut.appendValue(hStats[hIdx].displayLabel);
			alleleStatsOut.appendValue(hStats[hIdx].sampleCount);
			alleleStatsOut.appendValue(hStats[hIdx].maxReads);
			alleleStatsOut.appendValue(hStats[hIdx].maxReadsFraction);
		}
		alleleStatsOut.close();
	}

	
	private void outputSampleCalls (LocusResult locusResult, int targetIdx) throws AnalysisException {
		SampleLocusResult[] sResults = locusResult.sampleResults;
		
		Target t = locusResult.locus.getTargets()[targetIdx];
		File outTargetFolder = FileUtilities.checkFolder(targetsFolder, t.getName(), true);
		
		String[] colNames = {"Sample","Call","Alleles","AlleleReads"};
		TableOutput sampleCallsOut = new TableOutput (outTargetFolder, "CallsBySample-final.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < sResults.length; sIdx++) {
			SampleLocusResult sResult = sResults[sIdx];
			sampleCallsOut.newRow();
			sampleCallsOut.appendValue(sResult.sample.getName());
			SampleCall call = sResult.alignmentTargets[targetIdx].stringentCall;
			String callString = call.getCallString();
			boolean useLenientCall = call.isMissing();
			if (useLenientCall) {
				call = sResult.alignmentTargets[targetIdx].lenientCall;
				callString = "["+call.getCallString()+"]";
			}
			sampleCallsOut.appendValue(callString);
			sampleCallsOut.appendValue(call.aaAllele);
			sampleCallsOut.appendValue(call.aaAlleleSummary);
		}
		sampleCallsOut.close();
	}
	
	
	private void outputSampleCallsAllTargets (LocusResult[] locusResults) throws AnalysisException {
		// Get the file headers
		ArrayList<String> colNameList = new ArrayList<String>();
		colNameList.add("Sample");
		for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
			Target[] targets = locusResults[lIdx].locus.getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				colNameList.add(targets[tIdx].getName()+"["+targets[tIdx].getTargetAaSeq()+"]");
			}
		}
		String[] colNames = colNameList.toArray(new String[colNameList.size()]);

		TableOutput sampleCallsOut = new TableOutput (rootFolder, "CallsBySampleAllTargets.tab", colNames, bufferSize);
		TableOutput sampleNrefCallsOut = new TableOutput (rootFolder, "NrefCallsBySampleAllTargets.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			sampleCallsOut.newRow();
			sampleCallsOut.appendValue(samples[sIdx].getName());
			sampleNrefCallsOut.newRow();
			sampleNrefCallsOut.appendValue(samples[sIdx].getName());
			
			for (int lIdx = 0; lIdx < locusResults.length; lIdx++) {
				LocusResult locusResult = locusResults[lIdx];
				SampleLocusResult sResult = locusResult.sampleResults[sIdx];
				for (int tIdx = 0; tIdx < sResult.alignmentTargets.length; tIdx++) {
					SampleCall call = sResult.alignmentTargets[tIdx].stringentCall;
					String aaAllele = call.aaAllele;
					String aaNrefAllele = call.aaNrefAllele;
					boolean useLenientCall = call.isMissing();
					if (useLenientCall) {
						call = sResult.alignmentTargets[tIdx].lenientCall;
						if (!call.isMissing()) {
							aaAllele = call.aaAllele.toLowerCase();
							aaNrefAllele = call.aaNrefAllele.toLowerCase();							
						}
					}
					sampleCallsOut.appendValue(aaAllele);
					sampleNrefCallsOut.appendValue(aaNrefAllele);
				}
			}
		}
		sampleCallsOut.close();
		sampleNrefCallsOut.close();
	}



	/* ==========================================================
	 * Execution
	 * ==========================================================
	 */
	public static void main(String[] args) {
		if (args.length < 4) {
			System.out.println("Usage: org.cggh.bam.genotyping.reads.ReadsAlignmentAnalysis <configFile> <sampleListFile> <refFasta> <rootFolder>");
			return;
		}
		File configFile = new File(args[0]);
		File sampleListFile = new File(args[1]);
		File refFastaFile = new File(args[2]);
		File rootFolder = new File(args[3]);
		
		try {
			TargetsAlignmentAnalysis task = new TargetsAlignmentAnalysis(configFile, sampleListFile, refFastaFile, rootFolder);
			task.execute();
		} catch (Exception e) {
			log.error("Error: " + e);
			return;
		}
		log.info("Exiting");
	}

}
