package org.cggh.bam.codon;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.*;


public class SampleCodonAnalyzer {
	
	private Sample        sample;
	private TargetAnalysisConfig config;
	private TargetLocus[] loci;
	private SampleCaller  caller;

	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleCodonAnalyzer (TargetAnalysisConfig config, Sample sample) throws AnalysisException  {
		this.config = config;
		this.loci = config.getLoci();
		this.sample = sample;		
		this.caller = new SampleCaller(config);
		ReadsAlignment.configure(config);
	}

	public SampleResults analyzeSample () throws AnalysisException, IOException  {
		
		// Read the reads from the SAM file
		ReadsRetriever srr = new ReadsRetrieverFromAlignment (config);
		ArrayList<Read>[] mappedReadLists = srr.retrieveSampleReads(sample);
		
		// Analyze each locus
		SampleLocusResult[] locusResults = new SampleLocusResult[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			TargetLocus locus = loci[lIdx];
			Target[] targets = locus.getTargets();
			SampleTargetResult[] tResults = new SampleTargetResult[targets.length];

			// Handle the case where we have no mapped reads
			ArrayList<Read> mappedReadList = mappedReadLists[lIdx];
			if (mappedReadList.isEmpty()) {
				for (int tIdx = 0; tIdx < targets.length; tIdx++) {
					tResults[tIdx] = new SampleTargetResult(targets[tIdx], sample);
				}
				locusResults[lIdx] = new SampleLocusResult(locus, tResults, 0, 0);
				continue;
			}
			
			// Make an alignment and discard those reads that have too many differences from consensus
			Read[] sampleReads = mappedReadList.toArray(new Read[mappedReadList.size()]);
			ReadsAlignment ra = new ReadsAlignment(sample, locus, sampleReads);
			sampleReads = ra.getAlignedReads();
			
			// Genotype the targets for each read
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				AlignmentTarget target = (AlignmentTarget)targets[tIdx];
				TargetGenotyper tg = new TargetGenotyper (target, config);
				TargetGenotype[] targetGenos = tg.extractTargetGenotypes (sampleReads);
				
				LabelCounters ntAlleleCounters = new LabelCounters();
				int lowQualityCount = 0;
				for (int rIdx = 0; rIdx < targetGenos.length; rIdx++) {
					TargetGenotype geno = targetGenos[rIdx];
					if (geno.isValidGenotype()) {
						String ntAllele = geno.getNtGenotype();
						ntAlleleCounters.increment(ntAllele);
					} else if (geno instanceof TargetGenotype.LowQualityTargetGenotype) {
						lowQualityCount++;
					}
				}
				
				// Final step: make calls
				SampleCall sampleCall = caller.callSample(target, ntAlleleCounters);
				
				// Store target results
				tResults[tIdx] = new SampleTargetResult(target, sample, sampleCall, ntAlleleCounters, lowQualityCount);
			}
			// Store locus results
			locusResults[lIdx] = new SampleLocusResult(locus, tResults, sampleReads.length, ra.getMisalignedReads().length);
		}
		return new SampleResults(sample, locusResults);
	}
	

	/* ************************************************************
	 * Results container classes
	 * ************************************************************
	 */
	public static class SampleResults {
		private Sample              sample;
		private SampleLocusResult[] locusResults;
		
		public SampleResults(Sample sample, SampleLocusResult[] locusResults) {
			this.sample = sample;
			this.locusResults = locusResults;
		}
		
		public Sample getSample() {
			return sample;
		}

		public SampleLocusResult[] getLocusResults() {
			return locusResults;
		}
	}
	
	public static class SampleLocusResult {
		private TargetLocus          locus;
		private SampleTargetResult[] targetResults;
		private int                  alignedCount;
		private int                  misalignedCount;
		
		public SampleLocusResult(TargetLocus locus, SampleTargetResult[] targetResults, int alignedCount, int misalignedCount) {
			this.locus = locus;
			this.targetResults = targetResults;
			this.alignedCount = alignedCount;
			this.misalignedCount = misalignedCount;
		}

		public TargetLocus getLocus() {
			return locus;
		}

		public SampleTargetResult[] getTargetResults() {
			return targetResults;
		}
		
		public int getAlignedCount() {
			return alignedCount;
		}
		
		public int getMisalignedCount() {
			return misalignedCount;
		}
	}
		
}
