package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.sampleClass.AlleleClassGenotype.*;
import org.cggh.bam.sampleClass.AlleleClassTarget.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.TableOutput;

//import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class AlleleClassAnalyzer {
	
	//private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private SampleClassConfig  config;
	private Sample             sample;
	private AlleleClassLocus[] loci;
	private File outFolder;
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public AlleleClassAnalyzer (SampleClassConfig config, Sample sample, File outFolder) throws AnalysisException  {
		this.config = config;
		this.loci = (AlleleClassLocus[])config.getLoci();
		this.sample = sample;
		this.outFolder = outFolder;
	}

	public SampleAlleleClassResults analyzeSample () throws AnalysisException, IOException  {
		
		// Read the reads from the SAM file
		SampleReadsRetriever srr = new SampleReadsRetriever (config);
		ArrayList<MappedRead>[] mappedReadLists = srr.retrieveSampleReads(sample);
		//outputSampleReads (mappedReadLists);
		
		SampleAlleleClassLocusResult[] locusResults = new SampleAlleleClassLocusResult[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			AlleleClassLocus locus = loci[lIdx];
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			
			MappedRead[] sampleReads = mappedReadList.toArray(new MappedRead[mappedReadList.size()]);
			
			AlleleClassTarget[] targets = locus.getTargets();
			AlleleClassResult[] tResults = new AlleleClassResult[targets.length];
			
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
				tResults[tIdx] = new AlleleClassResult(target, alleleSetCounters, unlistedAlleleCounters.getSortedCounters());
			}
			locusResults[lIdx] = new SampleAlleleClassLocusResult(locus, tResults);
		}
		return new SampleAlleleClassResults(sample, locusResults);
	}
	
	protected void outputSampleReads (ArrayList<MappedRead>[] mappedReadLists) throws AnalysisException, IOException  {
		int idx = 0;
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			AlleleClassLocus locus = loci[lIdx];
			AlleleClassTarget[] targets = locus.getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				AlleleClassTarget target = targets[tIdx];
				ArrayList<MappedRead> reads = mappedReadLists[idx++];
				String[] headers = {"Sam"};
				TableOutput out = new TableOutput (outFolder, sample.getName()+"_"+locus.getName()+"_"+target.getName()+".reads.txt", headers, 64 * 1024);		
				for (int aIdx = 0; aIdx < reads.size(); aIdx++) {
					MappedRead r = reads.get(aIdx);
					out.newRow();
					out.appendValue(r.getSamString());
				}
				out.close();
			}
		}
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
		private AlleleClassResult[]   targetResults;
		
		public SampleAlleleClassLocusResult(AlleleClassLocus locus, AlleleClassResult[] targetResults) {
			this.locus = locus;
			this.targetResults = targetResults;
		}

		public AlleleClassLocus getLocus() {
			return locus;
		}

		public AlleleClassResult[] getTargetResults() {
			return targetResults;
		}
	}
}
