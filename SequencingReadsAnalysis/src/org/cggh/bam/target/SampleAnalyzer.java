package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.*;


public class SampleAnalyzer {
	
	private Sample        sample;
	private BamConfig     config;
	private TargetLocus[] loci;
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleAnalyzer (BamConfig config, Sample sample) throws AnalysisException  {
		this.config = config;		
		this.loci = config.getLoci();
		this.sample = sample;		
	}

	public SampleResults analyzeSample () throws AnalysisException, IOException  {
		
		// Read the reads from the SAM file
		SampleReadsRetriever srr = new SampleReadsRetriever (config);
		ArrayList<MappedRead>[] mappedReadLists = srr.retrieveSampleReads(sample);
		
		// Analyze each locus
		SampleLocusResult[] locusResults = new SampleLocusResult[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			TargetLocus locus = loci[lIdx];
			
			// Make an alignment and discard those reads that have too many differences from consensus
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			MappedRead[] sampleReads = mappedReadList.toArray(new MappedRead[mappedReadList.size()]);
			ReadsAlignment ra = new ReadsAlignment(sample, locus, sampleReads);
			sampleReads = ra.getAlignedReads();
			
			// Genotype the targets for each read
			Target[] targets = locus.getTargets();
			SampleTargetResult[] tResults = new SampleTargetResult[targets.length];
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				Target target = targets[tIdx];
				TargetGenotyper tg = new TargetGenotyper (target);
				TargetGenotype[] targetGenos = tg.extractTargetGenotypes (sampleReads);
				
				LabelCounters ntAlleleCounters = new LabelCounters();
				for (int rIdx = 0; rIdx < targetGenos.length; rIdx++) {
					TargetGenotype geno = targetGenos[rIdx];
					if (geno.isValidGenotype()) {
						String ntAllele = geno.getNtGenotype();
						ntAlleleCounters.increment(ntAllele);
					}
				}
				tResults[tIdx] = new SampleTargetResult(target, sample, ntAlleleCounters);
			}
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
