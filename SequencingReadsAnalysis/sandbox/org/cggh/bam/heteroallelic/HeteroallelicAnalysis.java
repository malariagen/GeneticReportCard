package org.cggh.bam.heteroallelic;

import org.cggh.bam.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;
import org.cggh.common.util.*;
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import org.apache.commons.logging.*;

public class HeteroallelicAnalysis extends SampleAnalysis {
	
	private static Log log = LogFactory.getLog((String)ClassUtilities.getCurrentClassName());
	
	public static final int MIN_PHRED_SCORE = 20;
	
	public static final int CALL_MISSING = 1;
	public static final int CALL_WT = 2;
	public static final int CALL_MUTANT = 3;
	public static final int CALL_HET = 4;
	
	private HeteroallelicConfig config;
	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

	private static final String[] CALL_FILE_HEADERS = new String[] { "Sample", "Locus", "Call", "Mutation", "MedianReadCount", "MeanReadCount" };
	private static final String[] MULTI_MUTANT_FILE_HEADERS = new String[] { "Sample", "Locus", "Alleles","ReadCount" };

	public HeteroallelicAnalysis(File configFile, File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException {
		super(refFastaFile, chrMapFile, outRootFolder);
		config = new HeteroallelicConfig(configFile);
	}

	public void analyzeSample(Sample sample) throws AnalysisException {
		log.info("Starting " + sample.getName());
		HeteroallelicLocus[] loci = config.getLoci();
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			try {
				analyzeSampleAtLocus(sample, loci[lIdx]);
			} catch (AnalysisException e) {
				log.info("Aborted analysis of " + sample.getName() + " at locus " + loci[lIdx].getName() + ": " + e);
				e.printStackTrace();
			}
		}
		log.info("Completed " + sample.getName());
	}

	public void analyzeSampleAtLocus(Sample sample, HeteroallelicLocus locus) throws AnalysisException {
		GenomeRegion locusregion = locus.getRegion();
		String chrName = locusregion.getChromosome();
		chrName = ChromosomeMap.getMappedChromosomeName((String) chrName, (String) sample.getBamChromosomeMap());
		Sequence chrRefSeq = ReferenceGenome.getChrSequence((String) chrName);
		int locusStartPos = locusregion.getStartPos();
		int locusEndPos = locusregion.getStopPos();
		String locusRefSeq = chrRefSeq.getData().substring(locusStartPos - 1, locusEndPos);
		if (locus.isReverse()) {
			locusRefSeq = SequenceUtilities.getReverseComplementSequence((String) locusRefSeq);
		}
		int locusSeqLen = locusRefSeq.length();
		int locusCodonCount = locusSeqLen / 3;
		char[] locusRefAminos = SequenceUtilities.translateNtSequence((String) locusRefSeq).toCharArray();

		AminoAlleleCounter[] codonCounters = new AminoAlleleCounter[locusCodonCount];
		for (int i = 0; i < locusCodonCount; i++) {
			codonCounters[i] = new AminoAlleleCounter();
		}
		ArrayList<MutantAllele> mutantAlleleList = new ArrayList<MutantAllele>();
		HashMap<String, MultipleMutant> multipleMutantTable = new HashMap<String, MultipleMutant>();
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		SAMRecordIterator it = samReader.query(chrName, locusStartPos, locusEndPos, false);
		while (it.hasNext()) {
			int trimRightIdx;
			int readRightOffset;
			int lastCodonIdx;
			int firstCodonIdx;
			int trimLeftIdx;
			int readEndPos;
			SAMRecord record = (SAMRecord) it.next();
			boolean isUngapped = (record.getCigarLength() == 1) && record.getCigarString().endsWith("M");
			if (!isUngapped) {
				continue;
			}
			int readLen = record.getReadLength();
			int readStartPos = record.getAlignmentStart();
			int readLeftOffset = readStartPos - locusStartPos;
			if (readLeftOffset < 0) {
				firstCodonIdx = 0;
				trimLeftIdx = -readLeftOffset;
				int trimmedLen = readLen - trimLeftIdx;
				if (trimmedLen < 3) {
					continue;
				}
			} else {
				firstCodonIdx = readLeftOffset / 3;
				int trimLeftLen = readLeftOffset % 3;
				trimLeftIdx = 0;
				if (trimLeftLen > 0) {
					trimLeftIdx += 3 - trimLeftLen;
					firstCodonIdx++;
				}
				int trimmedLen = readLen - trimLeftIdx;
				if (trimmedLen < 3)
					continue;
			}
			readEndPos = readStartPos + readLen - 1;
			if ((readRightOffset = locusEndPos - readEndPos) < 0) {
				lastCodonIdx = locusCodonCount - 1;
				trimRightIdx = readLen - 1 + readRightOffset;
				int trimmedLen = 1 + trimRightIdx - trimLeftIdx;
				if (trimmedLen < 3) {
					continue;
				}
			} else {
				lastCodonIdx = locusCodonCount - 1 - readRightOffset / 3;
				int trimRightLen = readRightOffset % 3;
				trimRightIdx = readLen - 1;
				if (trimRightLen > 0) {
					trimRightIdx -= 3 - readRightOffset % 3;
					--lastCodonIdx;
				}
				int trimmedLen = 1 + trimRightIdx - trimLeftIdx;
				if (trimmedLen < 3)
					continue;
			}
			String readSeq = record.getReadString();
			String readQ = record.getBaseQualityString();
			String trimmedReadSeq = readSeq.substring(trimLeftIdx, trimRightIdx + 1);
			String trimmedReadQ   = readQ.substring(trimLeftIdx, trimRightIdx + 1);
			if (locus.isReverse()) {
				trimmedReadSeq = SequenceUtilities.getReverseComplementSequence((String) trimmedReadSeq);
				trimmedReadQ = TextUtilities.reverse((String) trimmedReadQ);
				int firstCodonIdxTmp = locusCodonCount - lastCodonIdx - 1;
				lastCodonIdx = locusCodonCount - firstCodonIdx - 1;
				firstCodonIdx = firstCodonIdxTmp;
			}
			char[] readAminos = SequenceUtilities.translateNtSequence((String) trimmedReadSeq).toCharArray();
			boolean[] highQuality = getCodonQuality(trimmedReadQ);
			mutantAlleleList.clear();
			
			for (int i2 = 0; i2 < readAminos.length; i2++) {
				int codonIdx = firstCodonIdx + i2;
				if (codonIdx >= locusCodonCount)
					break;
				if (highQuality[i2]) {
					char amino = readAminos[i2];
					char refAmino = locusRefAminos[codonIdx];
					if (amino != refAmino) {
						mutantAlleleList.add(new MutantAllele(codonIdx, refAmino, amino));
					}
					codonCounters[codonIdx].increment(amino);
				}
			}
			if (mutantAlleleList.size() <= 1)
				continue;
			MultipleMutant mm = new MultipleMutant(mutantAlleleList);
			String mmLabel = mm.getLabel();
			MultipleMutant tableMm = (MultipleMutant) multipleMutantTable.get(mmLabel);
			if (tableMm == null) {
				multipleMutantTable.put(mmLabel, mm);
				tableMm = mm;
			}
			tableMm.readCounts++;
		}
		it.close();
		callSample(sample, locus, codonCounters, locusRefAminos);
		ArrayList<MultipleMutant> multipleMutantList = new ArrayList<MultipleMutant>();
		for (MultipleMutant mm : multipleMutantTable.values()) {
			if (mm.readCounts < 2)
				continue;
			multipleMutantList.add(mm);
		}
		if (!multipleMutantList.isEmpty()) {
			outputMultipleMutants(sample, locus, multipleMutantList);
		}
	}

	private boolean[] getCodonQuality(String qString) {
		int codonCount = qString.length() / 3;
		boolean[] hq = new boolean[codonCount];
		for (int i = 0; i < codonCount; i++) {
			int baseIdx = 3 * i;
			hq[i] =   (((qString.charAt(baseIdx)     - 33) >= 20)
					&& ((qString.charAt(baseIdx + 1) - 33) >= 20)
					&& ((qString.charAt(baseIdx + 2) - 33) >= 20));
		}
		return hq;
	}

	public void callSample(Sample sample, HeteroallelicLocus locus, AminoAlleleCounter[] codonCounters, char[] refAminos) throws AnalysisException {
		File outFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
		TableOutput alleleOut = new TableOutput(outFolder, String.valueOf(sample.getName()) + "." + locus.getName() + ".mutations.tab", 
				new String[] { "Sample", "Locus", "Codon", "Call", "Mutation", "TotalReadCount", "MutantReadCount", "MutantReadProp" }, 64 * 1024);
		int sampleCall = CALL_WT;
		String sampleMutation = null;
		int[] readCounts = new int[codonCounters.length];
		
		for (int cIdx = 0; cIdx < codonCounters.length; cIdx++) {
			CodonCall cc = callCodon(cIdx + locus.getStartCodon(), codonCounters[cIdx], refAminos[cIdx]);
			readCounts[cIdx] = cc.totalReads;
			
	    switch (cc.call) {
			case CALL_MISSING:
			case CALL_WT: {
				break;
			}
			case CALL_MUTANT: {
				switch (sampleCall) {
				case CALL_MISSING:
				case CALL_WT: {
					sampleCall = CALL_MUTANT;
					sampleMutation = cc.mutation;
					break;
				}
				case CALL_MUTANT:
				case CALL_HET:
					sampleCall = CALL_HET;
					sampleMutation = sampleMutation + "," + cc.mutation;
					break;
				}
				break;
			}
			case CALL_HET:
				sampleCall = CALL_HET;
				switch (sampleCall) {
				case CALL_MISSING:
				case CALL_WT:
					sampleMutation = cc.mutation;
					break;
				case CALL_MUTANT:
				case CALL_HET:
					sampleMutation = sampleMutation + "," + cc.mutation;
					break;
				}
				break;
			}
			if (cc.call == CALL_MUTANT || cc.call == CALL_HET) {
				alleleOut.newRow();
				alleleOut.appendValue(sample.getName());
				alleleOut.appendValue(locus.getName());
				alleleOut.appendValue(cIdx + locus.getStartCodon());
				alleleOut.appendValue(getCallString(cc.call));
				alleleOut.appendValue(cc.mutation);
				alleleOut.appendValue(cc.totalReads);
				alleleOut.appendValue(cc.mutReads);
				alleleOut.appendValue(cc.mutProp);
			}
		}
		alleleOut.close();
		Statistics stat = new Statistics(readCounts);
		TableOutput callOut = new TableOutput(outFolder, sample.getName() + "." + locus.getName() + ".calls.tab", CALL_FILE_HEADERS, 1024);
		callOut.newRow();
		callOut.appendValue(sample.getName());
		callOut.appendValue(locus.getName());
		callOut.appendValue(getCallString(sampleCall));
		callOut.appendValue(sampleMutation);
		callOut.appendValue(stat.getMedian());
		callOut.appendValue(stat.getMean());
		callOut.close();
	}

	private String getCallString(int call) {
		switch (call) {
		case CALL_MISSING: 
			return "MI";
		case CALL_WT:
			return "WT";
		case CALL_MUTANT:
			return "MU";
		case CALL_HET:
			return "HE";
		}
		return "-";
	}

	public CodonCall callCodon(int codonPos, AminoAlleleCounter codonCounter, char refAllele) throws AnalysisException {
		CodonCall cc = new CodonCall();
		AlleleCounter.AlleleCount[] counts = codonCounter.getSortedAlleleCounts();
		
		for (int i = 0; i < counts.length; i++) {
			int alleleNumReads = counts[i].getCount();
			if (alleleNumReads <= 2) {
				break;
			}
			char allele = counts[i].getAllele();
			if (allele == refAllele) {
				cc.isWt = true;
			} else {
				cc.mutReads += alleleNumReads;
			}
			cc.alleleCount++;
			cc.totalReads += alleleNumReads;
		}
		if (cc.totalReads >= 5) {
			if (cc.alleleCount == 1) {
				if (cc.isWt) {
					cc.call = CALL_WT;
				} else {
					cc.call = CALL_MUTANT;
					cc.mutation = "" + refAllele + codonPos + counts[0].getAllele();
					cc.mutReads = cc.totalReads;
					cc.mutProp = 1.0;
				}
			} else {
				cc.call = CALL_HET;
				char topAllele = counts[0].getAllele();
				char nonRefAllele = topAllele == refAllele ? counts[1].getAllele() : topAllele;
				cc.mutation = "" + refAllele + codonPos + nonRefAllele + '*';
				double wtReads = topAllele == refAllele ? counts[0].getCount() : counts[1].getCount();
				cc.mutProp = 1.0 - wtReads / (double) cc.totalReads;
			}
		}
		return cc;
	}

	public void outputMultipleMutants(Sample sample, HeteroallelicLocus locus, Collection<MultipleMutant> mutantsList) throws AnalysisException {
		File outFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
		TableOutput mmOut = new TableOutput(outFolder, String.valueOf(sample.getName()) + "." + locus.getName() + ".multipleMutants.tab", MULTI_MUTANT_FILE_HEADERS, 4096);
		for (MultipleMutant mm : mutantsList) {
			mmOut.newRow();
			mmOut.appendValue(sample.getName());
			mmOut.appendValue(locus.getName());
			mmOut.appendValue(mm.getLabel());
			mmOut.appendValue(mm.readCounts);
		}
		mmOut.close();
	}

	public void analyzeAllSampleResults(Sample[] samples) throws AnalysisException, IOException {
		HeteroallelicLocus[] loci = config.getLoci();
		
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			try {
				analyzeAllSamplesAtLocus(samples, loci[lIdx]);
			} catch (Exception e) {
				log.info("Aborted analysis of locus " + loci[lIdx].getName() + ": " + e);
				e.printStackTrace();
			}
		}
	}

	private void analyzeAllSamplesAtLocus(Sample[] samples, HeteroallelicLocus locus) throws AnalysisException, IOException {
		mergeResultFiles(samples, locus, ".calls", CALL_FILE_HEADERS);
		mergeResultFiles(samples, locus, ".multipleMutants", MULTI_MUTANT_FILE_HEADERS);
	}

	private void mergeResultFiles(Sample[] samples, HeteroallelicLocus locus, String filenameSuffix, String[] fieldHeaders) throws AnalysisException, IOException {
		TableOutput mergeOut = new TableOutput(outRootFolder, "AllSamples." + locus.getName() + filenameSuffix + ".tab", fieldHeaders, 1048576);
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
			File sampleFile = new File(sampleFolder, String.valueOf(sample.getName()) + "." + locus.getName() + filenameSuffix + ".tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				TableInput tif = new TableInput(sampleFile);
				tif.getFieldNames();
				try {
					String[] inFields;
					while ((inFields = tif.getNextValidLine()) != null) {
						mergeOut.newRow();
						for (int fIdx = 0; fIdx < fieldHeaders.length; fIdx++) {
							mergeOut.appendValue(inFields[fIdx + 1]);
						}
					}
				} finally {
					tif.close();
				}
			}
		}
		mergeOut.close();
	}

	public static class CodonCall {
		int call = CALL_MISSING;
		String mutation = "-";
		int mutReads = 0;
		int totalReads = 0;
		boolean isWt = false;
		int alleleCount = 0;
		double mutProp = 0.0;
	}

	public static class MultipleMutant {
		MutantAllele[] mutantAlleles;
		int readCounts;
		String label;

		public MultipleMutant(ArrayList<MutantAllele> mutantAlleleList) {
			mutantAlleles = mutantAlleleList.toArray(new MutantAllele[mutantAlleleList.size()]);
		}

		public String getLabel() {
			if (label == null) {
				StringBuffer sb = new StringBuffer();
				for (int i = 0; i < mutantAlleles.length; i++) {
					if (sb.length() > 0) {
						sb.append(',');
					}
					sb.append(mutantAlleles[i].toString());
				}
				label = sb.toString();
			}
			return label;
		}
	}

	public static class MutantAllele {
		int codonIdx;
		char refAllele;
		char allele;
		String label;
		int readCount;

		public MutantAllele(int codonIdx, char refAllele, char allele) {
			this.codonIdx = codonIdx;
			this.refAllele = refAllele;
			this.allele = allele;
		}

		public String toString() {
			if (label == null) {
				label = "" + refAllele + (codonIdx + 1) + allele;
			}
			return label;
		}
	}

	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		public static void main(String[] args) {
			if (args.length < 7) {
				log.error("Usage: org.cggh.bam.heteroallelic.HeteroallelicAnalysis$SingleSample <configFile> <sampleName> <bamFile> <chrMap> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);			log.info("ConfigFile: " + configFile.getAbsolutePath());
			String sampleId = args[1];						log.info("SampleId: " + sampleId);
			File sampleBamFile = new File(args[2]);			log.info("SampleBamFile: " + sampleBamFile.getAbsolutePath());
			String sampleChrMapName = args[3];				log.info("SampleChrMapName: " + sampleChrMapName);
			File refFastaFile = new File(args[4]);			log.info("RefFastaFile: " + refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[5]);			log.info("ChrMapFile: " + chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[6]);			log.info("RootFolder: " + rootFolder.getAbsolutePath());
			try {
				Sample sample = new Sample(sampleId, sampleBamFile, sampleChrMapName);
				HeteroallelicAnalysis task = new HeteroallelicAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
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
				log.error("Usage: org.cggh.bam.heteroallelic.HeteroallelicAnalysis$MultiSample <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);			log.info("ConfigFile: " + configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);		log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);			log.info("RefFastaFile: " + refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);			log.info("ChrMapFile: " + chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);			log.info("RootFolder: " + rootFolder.getAbsolutePath());
			
			
			int maxThreads = Integer.parseInt(System.getProperty("maxThreads", "0"));
			
			try {
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				HeteroallelicAnalysis task = new HeteroallelicAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
				multi.execute((SampleAnalysis) task);
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
				log.error("Usage: org.cggh.bam.heteroallelic.HeteroallelicAnalysis$MergeResults <configFile> <sampleListFile> <refFasta> <chrMapFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);			log.info("ConfigFile: " + configFile.getAbsolutePath());
			File sampleListFile = new File(args[1]);		log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File refFastaFile = new File(args[2]);			log.info("RefFastaFile: " + refFastaFile.getAbsolutePath());
			File chrMapFile = new File(args[3]);			log.info("ChrMapFile: " + chrMapFile.getAbsolutePath());
			File rootFolder = new File(args[4]);			log.info("RootFolder: " + rootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				HeteroallelicAnalysis task = new HeteroallelicAnalysis(configFile, refFastaFile, chrMapFile, rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
