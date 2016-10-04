package org.cggh.bam.target.alleleClasses;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.bam.target.alleleClasses.AlleleClassTarget.*;
import org.cggh.bam.target.alleleClasses.AlleleClassGenotype.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
//import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class SampleAlleleClassAnalyzer {
	
	//private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private Sample                sample;
	private AlleleClassLocus[] loci;
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleAlleleClassAnalyzer (AlleleClassLocus[] loci, Sample sample) throws AnalysisException  {
		this.loci = loci;
		this.sample = sample;		
	}

	public SampleAlleleClassResults analyzeSample () throws AnalysisException, IOException  {
		
		// **************************************
		// PART 1 - Reads analysis and Alignment
		// **************************************
		// Read the reads from the SAM file
		SampleReadsRetriever srr = new SampleReadsRetriever (loci);
		ArrayList<MappedRead>[] mappedReadLists = srr.retrieveSampleReads(sample);
		
		SampleAlleleClassLocusResult[] locusResults = new SampleAlleleClassLocusResult[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			AlleleClassLocus locus = loci[lIdx];
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			
			/*
			// Remove reads that do not cover the target
			MappedReadFilter filter1 = new MappedReadTargetQualityFilter(locus.targets);
			@SuppressWarnings("unused")
			ArrayList<MappedRead> discardedReadsList = filter1.filterReads(mappedReadList);
			*/
			MappedRead[] sampleReads = mappedReadList.toArray(new MappedRead[mappedReadList.size()]);
			
			AlleleClassTarget[] targets = locus.getTargets();
			SampleAlleleClassResult[] tResults = new SampleAlleleClassResult[targets.length];
			
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				AlleleClassTarget target = targets[tIdx];
				AlleleClassGenotyper tg = new AlleleClassGenotyper (target);
				TargetGenotype[] targetGenos = tg.extractTargetNtGenotypes (sampleReads);
				
				ClassAllele[] alleles = target.getAlleles();
				LabelCounter[] alleleSetCounters = new LabelCounter[alleles.length];
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					alleleSetCounters[aIdx] = new LabelCounter(alleles[aIdx].getName());
				}
				
				LabelCounters unlistedAlleleCounters = new LabelCounters();
				
				for (int rIdx = 0; rIdx < targetGenos.length; rIdx++) {
					TargetGenotype geno = targetGenos[rIdx];
					if (geno instanceof AlleleClassGenotype) {
						int alleleIdx = ((AlleleClassGenotype)geno).alleleIdx;
						alleleSetCounters[alleleIdx].increment();
					} else if (geno instanceof UnlistedAlleleTargetGenotype) {
						unlistedAlleleCounters.increment(geno.getNtGenotype());
					}
				}
				tResults[tIdx] = new SampleAlleleClassResult(target, alleleSetCounters, unlistedAlleleCounters.getSortedCounters());
			}
			locusResults[lIdx] = new SampleAlleleClassLocusResult(locus, tResults);
		}
		return new SampleAlleleClassResults(sample, locusResults);
	}
	

	/* ************************************************************
	 * Results container classes
	 * ************************************************************
	 */
	public static class SampleAlleleClassResults {
		private Sample              sample;
		private SampleAlleleClassLocusResult[] locusResults;
		
		public SampleAlleleClassResults(Sample sample, SampleAlleleClassLocusResult[] locusResults) {
			this.sample = sample;
			this.locusResults = locusResults;
		}
		
		public Sample getSample() {
			return sample;
		}

		public SampleAlleleClassLocusResult[] getLocusResults() {
			return locusResults;
		}
	}
	
	public static class SampleAlleleClassLocusResult {
		private AlleleClassLocus locus;
		private SampleAlleleClassResult[]   targetResults;
		
		public SampleAlleleClassLocusResult(AlleleClassLocus locus, SampleAlleleClassResult[] targetResults) {
			this.locus = locus;
			this.targetResults = targetResults;
		}

		public AlleleClassLocus getLocus() {
			return locus;
		}

		public SampleAlleleClassResult[] getTargetResults() {
			return targetResults;
		}
	}
}
