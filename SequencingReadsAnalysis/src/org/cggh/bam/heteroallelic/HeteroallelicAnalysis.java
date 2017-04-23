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
	
	public static final int MIN_CALL_COVERAGE = 5;
	public static final double MAX_READ_DIRECTION_PROPORTION = 0.85;
	
	public static final int MIN_PHRED_SCORE = 20;
	
	public static final int CALL_MISSING = 1;
	public static final int CALL_WT = 2;
	public static final int CALL_MUTANT = 3;
	public static final int CALL_HET = 4;
	
	private HeteroallelicConfig   config;
	private SamReaderFactory      samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	private SampleLocusLogFile    msgLog;
	private Genotyper genotyper = new Genotyper.GenotyperReadCountProportion(0.05); // 5% total reads is the min to call an allele
	

	private static final String[] CALL_FILE_HEADERS = new String[] { "Sample", "Locus", "Call", "Mutation", "MissingCodonCallsProp", "MedianReadCount", "MeanReadCount" };
	private static final String[] MULTI_MUTANT_FILE_HEADERS = new String[] { "Sample", "Locus", "Alleles", "ReadCount" };
	private static final String[] MESSAGE_LOG_FILE_HEADERS = new String[] { "Message", "Allele", "ReadCount" };

	public HeteroallelicAnalysis(File configFile, File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException {
		super(refFastaFile, chrMapFile, outRootFolder);
		config = new HeteroallelicConfig(configFile);
	}
	
	public void analyzeSample(Sample sample) throws AnalysisException {
		log.info("Starting " + sample.getName());
		HeteroallelicLocus[] loci = config.getLoci();
		File outFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
		
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			HeteroallelicLocus locus = loci[lIdx];
			msgLog = new SampleLocusLogFile(sample, locus, MESSAGE_LOG_FILE_HEADERS);
			
			try {
				// Analyze the sample in the locus region
				LocusRegion region = analyzeSampleLocusRegion(sample, locus);
				
				// Write out the mutant codon calls to file
				TableOutput alleleOut = new TableOutput(outFolder, sample.getName() + "." + locus.getName() + ".mutations.tab", 
						new String[] { "Sample", "Locus", "Codon", "Call", "Mutation", "TotalReadCount", "MutantReadCount", "MutantReadProp" }, 64 * 1024);
				int sampleCall = CALL_MISSING;
				String sampleMutation = null;
				int[] readCounts = new int[region.codonCount];
				
				int missingCount = 0;
				for (int cIdx = 0; cIdx < region.codonCount; cIdx++) {
					CodonCall cc = region.callCodon(cIdx);
					readCounts[cIdx] = cc.totalReads;
					
			        switch (cc.call) {
					case CALL_MISSING:
						missingCount++;
						break;
					case CALL_WT:
						if (sampleCall == CALL_MISSING) {
							sampleCall = CALL_WT;
							sampleMutation = "WT";
						}
						break;
					case CALL_MUTANT: 
						switch (sampleCall) {
						case CALL_MISSING:
						case CALL_WT:
							sampleCall = CALL_MUTANT;
							sampleMutation = cc.mutation;
							break;
						case CALL_MUTANT:
						case CALL_HET:
							sampleCall = CALL_HET;
							sampleMutation = (sampleMutation == null) ? cc.mutation : sampleMutation+","+cc.mutation;
							break;
						}
						break;
					case CALL_HET:
						sampleCall = CALL_HET;
						switch (sampleCall) {
						case CALL_MISSING:
						case CALL_WT:
							sampleMutation = cc.mutation;
							break;
						case CALL_MUTANT:
						case CALL_HET:
							sampleMutation = (sampleMutation == null) ? cc.mutation : sampleMutation+","+cc.mutation;
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
				
				double missingCallsProp = ((double)missingCount) / ((double)region.codonCount);

				// Write out the call for the sample to file
				Statistics stat = new Statistics(readCounts);
				TableOutput callOut = new TableOutput(outFolder, sample.getName() + "." + locus.getName() + ".calls.tab", CALL_FILE_HEADERS, 1024);
				callOut.setMaximumFractionDigits(2);
				callOut.newRow();
				callOut.appendValue(sample.getName());
				callOut.appendValue(locus.getName());
				callOut.appendValue(getCallString(sampleCall));
				callOut.appendValue(sampleMutation);
				callOut.appendValue(missingCallsProp);
				callOut.appendValue(stat.getMedian());
				callOut.appendValue(stat.getMean());
				callOut.close();
				
				// Report any multiple mutant reads observed in the sample
				MultipleMutant[] mutants = region.getMultipleMutants ();
				if (mutants.length > 0) {
					TableOutput mmOut = new TableOutput(outFolder, sample.getName() + "." + locus.getName() + ".multipleMutants.tab", MULTI_MUTANT_FILE_HEADERS, 4096);
					for (MultipleMutant mm : mutants) {
						mmOut.newRow();
						mmOut.appendValue(sample.getName());
						mmOut.appendValue(locus.getName());
						mmOut.appendValue(mm.getLabel());
						mmOut.appendValue(mm.readCounts);
					}
					mmOut.close();
				}
				
				// Report messages, if any were generated for the sample
				if (!msgLog.isEmpty()) {
					msgLog.saveFile(outFolder);
					msgLog.clear();
				}
				
			} catch (AnalysisException e) {
				log.info("Aborted analysis of " + sample.getName() + " at locus " + loci[lIdx].getName() + ": " + e);
				e.printStackTrace();
			}
		}
		
		log.info("Completed " + sample.getName());
	}
	
	public LocusRegion analyzeSampleLocusRegion(Sample sample, HeteroallelicLocus locus) throws AnalysisException {
		
		LocusRegion region = new LocusRegion (locus, sample);
		ArrayList<MutantAllele> mutantAlleleList = new ArrayList<MutantAllele>();
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		SAMRecordIterator it = samReader.query(region.chrName, region.startPos, region.endPos, false);
		while (it.hasNext()) {
			SAMRecord record = (SAMRecord) it.next();
			boolean isUngapped = (record.getCigarLength() == 1) && record.getCigarString().endsWith("M");
			if (!isUngapped) {
				continue;
			}
			int readLen = record.getReadLength();
			int readStartPos = record.getAlignmentStart();
			int readLeftOffset = readStartPos - region.startPos;
			int firstCodonIdx;
			int trimLeftIdx;
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
			
			int trimRightIdx;
			int readRightOffset;
			int lastCodonIdx;
			int readEndPos = readStartPos + readLen - 1;
			if ((readRightOffset = region.endPos - readEndPos) < 0) {
				lastCodonIdx = region.codonCount - 1;
				trimRightIdx = readLen - 1 + readRightOffset;
				int trimmedLen = 1 + trimRightIdx - trimLeftIdx;
				if (trimmedLen < 3) {
					continue;
				}
			} else {
				lastCodonIdx = region.codonCount - 1 - readRightOffset / 3;
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
				int firstCodonIdxTmp = region.codonCount - lastCodonIdx - 1;
				lastCodonIdx = region.codonCount - firstCodonIdx - 1;
				firstCodonIdx = firstCodonIdxTmp;
			}
			char[] readAminos = SequenceUtilities.translateNtSequence((String) trimmedReadSeq).toCharArray();
			boolean[] highQuality = getCodonQuality(trimmedReadQ);
			mutantAlleleList.clear();
			
			boolean isReverse = record.getReadNegativeStrandFlag();
			for (int i2 = 0; i2 < readAminos.length; i2++) {
				int codonIdx = firstCodonIdx + i2;
				if (codonIdx >= region.codonCount)
					break;
				if (highQuality[i2]) {
					char amino = readAminos[i2];
					char refAmino = region.refAminos[codonIdx];
					if (amino != refAmino) {
						mutantAlleleList.add(new MutantAllele(codonIdx+locus.getStartCodon(), refAmino, amino));
					}
					region.processReadAminoAllele (codonIdx, amino, isReverse);
				}
			}
			
			region.processReadForMultipleMutations (mutantAlleleList);
		}
		it.close();
		return region;
	}
	
	private class LocusRegion {
		HeteroallelicLocus locus;
		//Sample sample;
		String chrName;
		int    startPos;
		int    endPos;
		String ntRefSeq;
		char[] refAminos;
		int    codonCount;
		AminoAlleleCounter[] codonCounters;    // Counts the reads for each allele at each position
		AminoAlleleCounter[] revCodonCounters; // Counts the reverse-strand reads for each allele at each position
		HashMap<String, MultipleMutant> multipleMutantTable;
		
		public LocusRegion (HeteroallelicLocus locus, Sample sample) throws AnalysisException {
			this.locus = locus;
			// Get the locus ref NT sequence
			GenomeRegion locusregion = locus.getRegion();
			chrName = locusregion.getChromosome();
			chrName = ChromosomeMap.getMappedChromosomeName(chrName, sample.getBamChromosomeMap());
			Sequence chrRefSeq = ReferenceGenome.getChrSequence(chrName);
			startPos = locusregion.getStartPos();
			endPos = locusregion.getStopPos();
			ntRefSeq = chrRefSeq.getData().substring(startPos - 1, endPos);
			if (locus.isReverse()) {
				ntRefSeq = SequenceUtilities.getReverseComplementSequence(ntRefSeq);
			}
			
			// Get the ref aminos
			refAminos = SequenceUtilities.translateNtSequence(ntRefSeq).toCharArray();
			codonCount = refAminos.length;

			codonCounters = new AminoAlleleCounter[codonCount];
			revCodonCounters = new AminoAlleleCounter[codonCount]; // Counts 
			for (int i = 0; i < codonCount; i++) {
				codonCounters[i] = new AminoAlleleCounter();
				revCodonCounters[i] = new AminoAlleleCounter();
			}
			multipleMutantTable = new HashMap<String, MultipleMutant>();
		}
		
		
		public void processReadAminoAllele (int codonIdx, char amino, boolean isReverse) {
			codonCounters[codonIdx].increment(amino);
			if (isReverse) {
				revCodonCounters[codonIdx].increment(amino);
			}
		}
		
		public void processReadForMultipleMutations (ArrayList<MutantAllele> mutantAlleleList) {
			if (mutantAlleleList.size() <= 1)
				return;
			MultipleMutant mm = new MultipleMutant(mutantAlleleList);
			String mmLabel = mm.getLabel();
			MultipleMutant tableMm = (MultipleMutant) multipleMutantTable.get(mmLabel);
			if (tableMm == null) {
				multipleMutantTable.put(mmLabel, mm);
				tableMm = mm;
			}
			tableMm.readCounts++;
		}
		
		public MultipleMutant[] getMultipleMutants () {
			ArrayList<MultipleMutant> multipleMutantList = new ArrayList<MultipleMutant>();
			for (MultipleMutant mm : multipleMutantTable.values()) {
				if (mm.readCounts < 2) {  // Remove singletons
					continue;
				}
				multipleMutantList.add(mm);
			}
			MultipleMutant[] result = new MultipleMutant[multipleMutantList.size()];
			return multipleMutantList.toArray(result);
		}
		
		public CodonCall callCodon(int cIdx) throws AnalysisException {
			
			int codonPos = cIdx + locus.getStartCodon();
			AminoAlleleCounter codonCounter = codonCounters[cIdx];
			AminoAlleleCounter revCodonCounter = revCodonCounters[cIdx];
			char refAllele = refAminos[cIdx];
			CodonCall cc = new CodonCall(codonPos, codonCounter, revCodonCounter, refAllele);
			return cc;
		}
	}


	public class CodonCall {

		int    call = CALL_MISSING;
		
		int    totalReads = 0;
		int    alleleCount = 0;
		
		String mutation = "-";
		int    mutReads = 0;
		double mutProp = 0.0;
		
		public CodonCall(int codonPos, AminoAlleleCounter codonCounter, AminoAlleleCounter revCodonCounter, char refAllele) {

			AlleleCounter.AlleleCount[] counts = codonCounter.getSortedAlleleCounts();
			
			// Fewer than 5 reads: missing
			totalReads = codonCounter.getCumulativeCount();
			if (totalReads < MIN_CALL_COVERAGE) {
				return;
			}
			
			// Count the alleles
			for (int i = 0; i < counts.length; i++) {
				int alleleReads = counts[i].getCount();
				boolean isValidAllele = genotyper.isValidAllele(alleleReads, totalReads);
				if (!isValidAllele) {
					break;
				}
				alleleCount++;
			}
			
			// One allele: it's either WT or mutant
			if (alleleCount == 1) {
				char allele = counts[0].getAllele();
				if (allele == refAllele) {
					call = CALL_WT;
				} else {
					call = CALL_MUTANT;
					mutation = "" + refAllele + codonPos + allele;
					mutReads = totalReads;
					mutProp = 1.0;
				}
				return;
			}
				
			// Multiple alleles: probably a het call, but we have to disregard errors and artefacts
			char allele0 = counts[0].getAllele();
			char allele1 = counts[1].getAllele();
			if ((allele0 != refAllele) && (allele1 != refAllele)) {
				// Non-biallelic call - give up on this one, at least for now - return missing
				msgLog.addMessage(new String[] { "MULTIALLELIC_CODON_WITHOUT_REF", ""+refAllele+codonPos+allele0+','+refAllele+codonPos+allele1, Integer.toString(totalReads)});
				return;
			}
			
			int nrefIdx = (allele0 == refAllele) ? 1 : 0;
			char nrefAllele = counts[nrefIdx].getAllele();
			
			totalReads = counts[0].getCount() + counts[1].getCount();
			int nrefReads = counts[nrefIdx].getCount();
			
			int refRevReads = revCodonCounter.getCount(refAllele);
			int nrefRevReads = revCodonCounter.getCount(nrefAllele);
			int totalRevReads = refRevReads + nrefRevReads;
			
			// Test if the alternative allele is the product of many reads on the same strad (likely artefact)
			boolean validNref = isValidNonRefAllele (nrefReads, nrefRevReads, totalReads, totalRevReads);
			if (!validNref) {
				call = CALL_WT;
				msgLog.addMessage(new String[] { "NONREF_ALLELE_FROM_SINGLE_STRAND_IN_HET", ""+refAllele+codonPos+nrefAllele+'*', Integer.toString(nrefReads) });
				return;
			}
			
			// Proper het
			call = CALL_HET;
			mutation = "" + refAllele + codonPos + nrefAllele + '*';
			mutReads = nrefReads;
			mutProp = ((double) nrefReads) / ((double) totalReads);
		}

		boolean isValidNonRefAllele (int nrefReads, int nrefRevReads, int totalReads, int totalRevReads) {
			// Discard if >85% of reads are in one direction 
			double directionProp = ((double) nrefRevReads) / ((double) nrefReads);
			return ((directionProp <= MAX_READ_DIRECTION_PROPORTION) 
			     && (directionProp >= (1.0 - MAX_READ_DIRECTION_PROPORTION)));
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

	/* ==========================================================
	 * Result and Utility classes
	 * ==========================================================
	 */
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
		int codonPos;
		char refAllele;
		char allele;
		String label;
		int readCount;

		public MutantAllele(int codonPos, char refAllele, char allele) {
			this.codonPos = codonPos;
			this.refAllele = refAllele;
			this.allele = allele;
		}

		public String toString() {
			if (label == null) {
				label = "" + refAllele + codonPos + allele;
			}
			return label;
		}
	}

	private static final String[] SAMPLE_LOCUS_HEADERS = new String[] { "Sample", "Locus"};
	private class SampleLocusLogFile extends MessageLogFile {
		String[] sampleLocusNames;
		
		public SampleLocusLogFile(Sample sample, HeteroallelicLocus locus, String[] headers) {
			super (TextUtilities.mergeStringLists(SAMPLE_LOCUS_HEADERS, headers));
			sampleLocusNames = new String[] { sample.getName(), locus.getName() };
		}

		@Override
		public void addMessage (String[] msgFields) {
			super.addMessage(TextUtilities.mergeStringLists(sampleLocusNames, msgFields));
		}
		
		public void saveFile (File outFolder) throws AnalysisException {
			String filename = sampleLocusNames[0] + "." + sampleLocusNames[1] + ".messages.tab";
			saveFile(outFolder, filename);
		}
	}

	/* ==========================================================
	 * Result Merging
	 * ==========================================================
	 */
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
		mergeResultFiles(samples, locus, ".calls");
		mergeResultFiles(samples, locus, ".mutations");
		mergeResultFiles(samples, locus, ".multipleMutants");
		mergeResultFiles(samples, locus, ".messages");
	}

	private void mergeResultFiles(Sample[] samples, HeteroallelicLocus locus, String filenameSuffix) throws AnalysisException, IOException {
		//TableOutput mergeOut = new TableOutput(outRootFolder, "AllSamples." + locus.getName() + filenameSuffix + ".tab", fieldHeaders, 1048576);
		TableOutput mergeOut = null;
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder(outRootFolder, sample.getName(), true);
			File sampleFile = new File(sampleFolder, String.valueOf(sample.getName()) + "." + locus.getName() + filenameSuffix + ".tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
			} else {
				TableInput tif = new TableInput(sampleFile);
				// Get the headers, removing the "Num" automatic field at the start
				String[] fieldHeaders = tif.getFieldNames();
				fieldHeaders = Arrays.copyOfRange(fieldHeaders, 1, fieldHeaders.length); 
				if (mergeOut == null) {
					mergeOut = new TableOutput(outRootFolder, "AllSamples." + locus.getName() + filenameSuffix + ".tab", fieldHeaders, 1048576);
				}
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
		if (mergeOut != null) {
			mergeOut.close();
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
