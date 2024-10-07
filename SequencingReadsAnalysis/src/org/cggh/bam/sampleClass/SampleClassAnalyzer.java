package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.genotyping.*;
import org.cggh.bam.sampleClass.ClassTarget.*;
import org.cggh.bam.sampleClass.ClassTargetGenotyper.*;
import org.cggh.bam.target.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;

import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class SampleClassAnalyzer {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private SampleClassConfig config;
	private Genotyper         genotyper;
	private ClassLocus[]      loci;
	private ClassTarget[]     allTargets;
	private File              outFolder;
	private double            minTargetCallProp;
	private double            minSpecificTargetCallProp;

	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleClassAnalyzer (SampleClassConfig config, File outFolder) throws AnalysisException  {
		this.config = config;
		this.outFolder = outFolder;
		this.genotyper = new Genotyper(config);
		this.minTargetCallProp = config.getMinCalledTargetProp();
		this.minSpecificTargetCallProp = config.getMinCalledSpecificTargetProp();

		this.loci = (ClassLocus[])config.getLoci();
		ArrayList<ClassTarget> tList = new ArrayList<ClassTarget>();
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			ClassTarget[] targets = loci[lIdx].getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				tList.add(targets[tIdx]);
			}
		}
		this.allTargets = tList.toArray(new ClassTarget[tList.size()]);
	}

	public SampleCall analyzeSample (Sample sample) throws AnalysisException, IOException  {
		
		
		// Read the reads from the SAM file
		ReadsRetriever srr = new ReadsRetrieverFromAlignment (config);
		ArrayList<Read>[] mappedReadLists = srr.retrieveSampleReads(sample);
		//outputSampleReads (sample, mappedReadLists);
		
		int tarIdx = 0;
		SampleTargetResult[] targetResults = new SampleTargetResult[allTargets.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			// Get all the reads for this locus
			ArrayList<Read> mappedReadList = mappedReadLists[lIdx];
			Read[] sampleReads = mappedReadList.toArray(new Read[mappedReadList.size()]);
			
			// Get the targets to be typed
			ClassTarget[] targets = loci[lIdx].getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				ClassTarget target = targets[tIdx];
				
				// Determine the target genotypes, one per read. These will be objects of a subclass of TargetGenotype.
				// The ones to be counted are of type ClassTargetGenotype; the others will be errors, low quality or unknown sequences.
				ClassTargetGenotyper tg = new ClassTargetGenotyper (target, config);
				TargetGenotype[] targetGenos = tg.extractTargetNtGenotypes (sampleReads);
				//outputTargetReadGenos (sample, loci[lIdx], target, targetGenos, sampleReads);
				
				// The ClassAllele is the class-named set of sequences that are matched for identifying the class
				// e.g. Pf@TATT|TAGT is one ClassAllele, where the class label is Pf and the sequences TATT and TAGT.
				// The alleleCounters will count the reads for each ClassAllele
				ClassAllele[] alleles = target.getAlleles();
				LabelCounter[] alleleCounters = new LabelCounter[alleles.length];
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					alleleCounters[aIdx] = new LabelCounter(alleles[aIdx].getName());
				}
				
				// Count the number of reads for each class allele (i.e. the targetGenos of class ClassTargetGenotype), 
				// and put unlisted alleles in separate counters (to inform analysis if there are other noteworthy alleles not currently matched)
				LabelCounters unlistedAlleleCounters = new LabelCounters();
				for (int rIdx = 0; rIdx < targetGenos.length; rIdx++) {
					TargetGenotype geno = targetGenos[rIdx];
					if (geno instanceof ClassTargetGenotype) {
						int alleleIdx = ((ClassTargetGenotype)geno).alleleIdx;
						alleleCounters[alleleIdx].increment();
					} else if (geno instanceof UnlistedAlleleTargetGenotype) {
						unlistedAlleleCounters.increment(geno.getNtGenotype());
					}
				}
				
				// We have all the read counts, now determine which ClassAllele(s) could be called at this target (discarding potential errors)
				String[] classTargetCalls = new String[alleles.length];
				int totalReads = 0;
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					totalReads += alleleCounters[aIdx].getCount();
				}
				boolean isValidCall = genotyper.hasSufficientReads(totalReads);
				
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					if (isValidCall) {
						boolean isValidAllele = genotyper.isValidAllele(alleleCounters[aIdx].getCount(), totalReads);
						classTargetCalls[aIdx] = isValidAllele ? alleles[aIdx].getName() : "-";
					} else {
						classTargetCalls[aIdx] = "-";
					}
				}
				targetResults[tarIdx++] = new SampleTargetResult(target, classTargetCalls, alleleCounters, unlistedAlleleCounters.getSortedCounters());
			}
		}
		
		// Make a call based on the target results
		// Get the list of sample classes we're trying to assign, and index them with a lookup table
		String[] classes = config.getClasses();
		HashMap<String,Integer> classesIdxTable = new HashMap<String,Integer>();
		for (int cIdx = 0; cIdx < classes.length; cIdx++) {
			classesIdxTable.put(classes[cIdx], cIdx);
		}
		
		// Process one class at a time to get the overall call for the sample
		int[] specificAlleleCounts     = new int[classes.length];
		int[] promiscuousAllelesCounts = new int[classes.length];
		for (int tIdx = 0; tIdx < allTargets.length; tIdx++) {
			ClassTarget target = allTargets[tIdx];
			SampleTargetResult tResult = targetResults[tIdx];
			
			String[] targetCalls = tResult.getClassTargetCalls();
			for (int tcIdx = 0; tcIdx < targetCalls.length; tcIdx++) {
				String alleleLabel = targetCalls[tcIdx];
				if ("-".equals(alleleLabel)) {
					continue;
				}
				String[] alleleClasses = alleleLabel.split("\\|");
				boolean isPromiscuous = alleleClasses.length > 1;
				for (String alleleClass : alleleClasses) {
					Integer idxObj = classesIdxTable.get(alleleClass);
					if (idxObj == null) {
						log.error("Allele '"+alleleClass+"' found at target "+target.getName()+" in sample "+ sample.getName()+ " is not a valid class");
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

		// Do the calling.
		// Check there are enough targets called, and enough targets where the class is called specifically (i.e. not a promiscuous target allele)
		String call = null;
		double minTargetCall         = minTargetCallProp * (double)allTargets.length;
		double minSpecificTargetCall = minSpecificTargetCallProp * (double)allTargets.length;
		for (int cIdx = 0; cIdx < classes.length; cIdx++) {
			boolean hasEnoughSpecificCalls = specificAlleleCounts[cIdx] >= minSpecificTargetCall;
			if (hasEnoughSpecificCalls) {
				int targetCallCount = specificAlleleCounts[cIdx] + promiscuousAllelesCounts[cIdx];
				boolean hasEnoughCalls = targetCallCount >= minTargetCall;
				if (hasEnoughCalls) {
					call = (call == null) ? classes[cIdx] : call+","+classes[cIdx];
				}
			}
		}
		if (call == null) {
			call = "-";
		}
		return new SampleCall(sample, call, targetResults);
	}

	
	protected void outputSampleReads (Sample sample, ArrayList<Read>[] mappedReadLists) throws AnalysisException, IOException  {
		int idx = 0;
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			ClassLocus locus = loci[lIdx];
			ClassTarget[] targets = locus.getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				ClassTarget target = targets[tIdx];
				ArrayList<Read> reads = mappedReadLists[idx++];
				String[] headers = {"Sam"};
				TableOutput out = new TableOutput (outFolder, sample.getName()+"_"+locus.getName()+"_"+target.getName()+".reads.txt", headers, 64 * 1024);		
				for (int aIdx = 0; aIdx < reads.size(); aIdx++) {
					Read r = reads.get(aIdx);
					out.newRow();
					out.appendValue(r.getSamString());
				}
				out.close();
			}
		}
	}
	
	protected void outputTargetReadGenos (Sample sample, Locus locus, ClassTarget target, TargetGenotype[] targetGenos, Read[] sampleReads) throws AnalysisException, IOException  {
		String[] headers = {"Read", "Genotype"};
		TableOutput out = new TableOutput (outFolder, sample.getName()+"_"+locus.getName()+"_"+target.getName()+".readsGenos.txt", headers, 64 * 1024);	
		
		for (int aIdx = 0; aIdx < sampleReads.length; aIdx++) {
			out.newRow();
			out.appendValue(sampleReads[aIdx].getId());
			out.appendValue(targetGenos[aIdx].getName());
		}
		out.close();
	}
	
}
