package org.cggh.bam.target.amino;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.*;


public class AminoSampleAnalyzer {
	
	private Sample             sample;
	private AminoTargetLocus[] loci;
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public AminoSampleAnalyzer (AminoTargetLocus[] loci, Sample sample) throws AnalysisException  {
		this.loci = loci;
		this.sample = sample;		
	}

	public SampleResults analyzeSample () throws AnalysisException, IOException  {
		
		// **************************************
		// PART 1 - Reads analysis and Alignment
		// **************************************
		// Read the reads from the SAM file
		SampleReadsRetriever srr = new SampleReadsRetriever (loci);
		ArrayList<MappedRead>[] mappedReadLists = srr.retrieveSampleReads(sample);
		
		SampleLocusResult[] locusResults = new SampleLocusResult[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			AminoTargetLocus locus = loci[lIdx];
			
			// Make an alignment and discard those reads that have too many differences from consensus
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			MappedRead[] sampleReads = mappedReadList.toArray(new MappedRead[mappedReadList.size()]);
			
			// Remove reads that do not align well with the rest of the reads at this locus
			ReadsAlignment ra = new ReadsAlignment(sample, locus, sampleReads);
			sampleReads = ra.getAlignedReads();
			
			// Genotype the targets for each read
			AminoTarget[] targets = locus.getTargets();
			SampleTargetResult[] tResults = new SampleTargetResult[targets.length];
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				AminoTarget target = targets[tIdx];
				TargetGenotyper tg = new AminoTargetGenotyper (target);
				TargetGenotype[] targetGenos = tg.extractTargetGenotypes (sampleReads);
				
				LabelCounters ntAlleleCounters = new LabelCounters();
				for (int rIdx = 0; rIdx < targetGenos.length; rIdx++) {
					TargetGenotype geno = targetGenos[rIdx];
					if (geno instanceof AminoTargetGenotype) {
						String ntAllele = geno.getNtGenotype();
						ntAlleleCounters.increment(ntAllele);
					}
				}
				tResults[tIdx] = new SampleTargetResult(target, ntAlleleCounters);
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
		private AminoTargetLocus     locus;
		private SampleTargetResult[] targetResults;
		private int                  alignedCount;
		private int                  misalignedCount;
		
		public SampleLocusResult(AminoTargetLocus locus, SampleTargetResult[] targetResults, int alignedCount, int misalignedCount) {
			this.locus = locus;
			this.targetResults = targetResults;
			this.alignedCount = alignedCount;
			this.misalignedCount = misalignedCount;
		}

		public AminoTargetLocus getLocus() {
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
		
	public static class SampleTargetResult {
		private AminoTarget    target;
		private LabelCounters  ntAlleleCounters;
		//private SampleCall     stringentCall;
		//private SampleCall     lenientCall;
		
		public SampleTargetResult(AminoTarget target, LabelCounters ntCounters) {
			this.target = target;
			this.ntAlleleCounters = ntCounters;
		}

		public AminoTarget getTarget() {
			return target;
		}

		public LabelCounters getNtAlleleCounters() {
			return ntAlleleCounters;
		}
		
		public void filterAlleleCounters (LabelCounterFilter filter) {
			ntAlleleCounters.filterCounters(filter);
		}
		
		public void cleanupTargetAlleles () {
			// Clean up the target sequence
			filterAlleleCounters (new LabelCounterFilterByLabelSubstring("N"));  // Remove sequences with undetermined nucleotide
			filterAlleleCounters (new LabelCounterFilterByMinCount(2));          // Remove Singletons			
		}
		/*
		public void callGenotypes () {
			stringentCall = new SampleCall (ntAlleleCounters, target, 5, 2);
			lenientCall = new SampleCall (ntAlleleCounters, target, 2, 2);
		}
		*/
		
	}
}
