package org.cggh.bam.barcode;

import org.apache.commons.logging.*;
import org.cggh.bam.*;
import org.cggh.common.counters.*;
import org.cggh.common.counters.AlleleCounter.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.*;
import org.cggh.common.textStore.*;

import htsjdk.tribble.readers.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;

import java.io.*;
import java.util.*;


public class BarcodeFromVcfAnalysis extends SampleAnalysis {
	
	public static final int MIN_READS_FOR_GENOTYPE = 5;

	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private GenomePosition[] snps;
	private HashMap<String,Integer> snpMap;

	public BarcodeFromVcfAnalysis (File outRootFolder, File snpListFile) throws AnalysisException {
		super (outRootFolder);
		snps = loadPositions (snpListFile);
		snpMap = new HashMap<String,Integer>(snps.length);
		for (int i = 0; i < snps.length; i++) {
			snpMap.put(snps[i].getName(), i);
		}
	}
	
	
	private Genotyper gt = new Genotyper.GenotyperReadCountProportion(0.05);
	
	public void analyzeSample(Sample sample) throws AnalysisException {
		String sampleName = sample.getName();
		File sampleFolder = getSampleSubfolder(outRootFolder, sampleName, false);
		File vcfFile = new File (sampleFolder, sampleName+".vcf.gz");
		
		Genotype[] genos = new Genotype[snps.length];
		NtAlleleCounter counter = new NtAlleleCounter();
		
		// Create the VCF reader
		InputStream vcfInStream = VcfFileUtils.createVcfFileReader (vcfFile);
		LineReader r = LineReaderUtil.fromBufferedStream(vcfInStream);
		LineIteratorImpl lineIt = new LineIteratorImpl(r);
		VCFCodec codec = new VCFCodec();
		
		// Read the sample IDs
		VCFHeader header = (VCFHeader)codec.readActualHeader(lineIt);
		List<String> sampleNamesList = header.getGenotypeSamples();
		String vcfSampleName = sampleNamesList.get(0);
		if (!vcfSampleName.equals(sampleName)) {
			throw new AnalysisException("Sample name in file " + vcfFile.getAbsolutePath() + " does not match: found '"+ vcfSampleName + "' instead of '" + sampleName + "'");
		}

		while (lineIt.hasNext()) {
			VariantContext ctx = codec.decode(lineIt.next());
			String chr = ctx.getContig();
			int pos = ctx.getStart();
			String snpName = chr + ':' + pos;
			Integer snpIdxObj = snpMap.get(snpName);
			if (snpIdxObj == null) {
				continue;
			}
			int snpIdx = snpIdxObj.intValue();
						
			// Read the read counts for all alleles
			counter.reset();
			int[] readCounts = ctx.getGenotype(0).getAD();

			// Get all alleles
			byte[] refBases = ctx.getReference().getBases();
			char refChar = Character.toUpperCase((char)refBases[0]);
			
			List<Allele> alleleList = ctx.getAlleles();
			for (int i = 0; i < alleleList.size(); i++) {
				Allele allele = alleleList.get(i);
				if (allele.isSymbolic()) {
					continue;  // Ignore
				}
				
				byte[] bases = allele.getBases();
				char alleleChar = Character.toUpperCase((char)bases[0]);
				counter.setCount(alleleChar, readCounts[i]);
			}
			
			genos[snpIdx] = new Genotype(counter, refChar);
		}
		r.close();
		
		// Write out the results
		outputSampleResults (sampleName, genos);
	}
		
	private class Genotype {

		private int    call = 0;
		private char   genotype = 'X';
		private char   topAllele = '-';
		private double nraf = Double.NaN;
		private int    refCount = 0;
		private int    nrefCount = 0;
		private String alleleCountStr = "-";
		

		public Genotype (NtAlleleCounter counter, char ref) {
			AlleleCount[] aCounts = counter.getSortedAlleleCounts();
			int totalReads = counter.getCumulativeCount();
			if (totalReads >= MIN_READS_FOR_GENOTYPE) {
				for (int i = 0; i < aCounts.length; i++) {
					AlleleCount c = aCounts[i];
					if (!gt.isValidAllele(c.getCount(), totalReads)) {
						break;
					}
					boolean isRef = (c.getAllele() == ref);
					if (i == 0) {
						topAllele = c.getAllele();
						call = isRef ? 1 : 2;
					} else {
						call = 3;
					}
					if (isRef) {
						this.refCount = c.getCount();
					} else {
						if (this.nrefCount == 0) {  // Only use one nref allele
							this.nrefCount = c.getCount();
						}
					}
				}				
			}
			
			
			switch (call) {
			case 1: 
				nraf = 0.0;
				genotype = topAllele;
				break;
			case 2: 
				nraf = 1.0;
				genotype = topAllele;
				break;
			case 3:
				nraf = ((double)nrefCount)/((double)(nrefCount+refCount));
				genotype = 'N';
				break;
			}
			alleleCountStr = AlleleCounter.makeAlleleCountString (aCounts);
		}

		protected Genotype(int call, char topAllele, int refCount, int nrefCount, String alleleCountStr) {
			this.call = call;
			this.topAllele = topAllele;
			this.refCount = refCount;
			this.nrefCount = nrefCount;
			this.alleleCountStr = alleleCountStr;
			
			switch (call) {
			case 0:
				nraf = Double.NaN;
				genotype = 'X';
				break;
			case 1: 
				nraf = 0.0;
				genotype = topAllele;
				break;
			case 2: 
				nraf = 1.0;
				genotype = topAllele;
				break;
			case 3:
				nraf = ((double)nrefCount)/((double)(nrefCount+refCount));
				genotype = 'N';
				break;
			}
		}
	}

	
	
	private void outputSampleResults(String sampleName, Genotype[] genos) throws AnalysisException {
		File sampleFolder = getSampleSubfolder(outRootFolder, sampleName, false);
	    String[] colNames = new String[] {"Chr","Pos","Genotype","TopAllele","Call","Nraf","Ref","Nref","Counts"};
		OutputTextStore outStore = new UncompressedOutputTextStore(sampleFolder, sampleName+".genos.tab");
		TableOutput out = new TableOutput (outStore, colNames, 1024 * 1024);
		
		for (int i = 0; i < genos.length; i++) {
			out.newRow();
			GenomePosition snp = snps[i];
			out.appendValue(snp.getChromosome());
			out.appendValue(snp.getPos());
			Genotype geno = genos[i];
			if (geno == null) {
				out.appendValue('X');
				out.appendBlankValue();
				out.appendValue(0);
				out.appendBlankValue();
				out.appendBlankValue();
				out.appendBlankValue();
				out.appendBlankValue();
			} else {
				out.appendValue(geno.genotype);
				out.appendValue(geno.topAllele);
				out.appendValue(geno.call);
				out.appendValue(geno.nraf);
				out.appendValue(geno.refCount);
				out.appendValue(geno.nrefCount);
				out.appendValue(geno.alleleCountStr);
			}
		}
		out.close();
	}
	

	/* ==========================================================
	 * Aggregation
	 * ==========================================================
	 */
	public void analyzeAllSampleResults (Sample[] samples) throws AnalysisException {
		Genotype[][] allGenos = new Genotype[samples.length][];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			String sampleName = samples[sIdx].getName();
			File sampleFolder = getSampleSubfolder(outRootFolder, sampleName, false);
			if (sampleFolder == null) {
				continue;
			}
			File resultFile = new File (sampleFolder, sampleName+".genos.tab");
			if (!resultFile.canRead()) {
				continue;
			}
			allGenos[sIdx] = loadSampleResults (resultFile);
		}
		
		// Output aggregated results
		outputBarcodes (samples, allGenos);
		outputAlleleCounts (samples, allGenos);
		outputNrafs (samples, allGenos);
	}
	
	private void outputBarcodes (Sample[] samples, Genotype[][] allGenos) throws AnalysisException {
		String[] colNames = new String[] {"Sample","MissingCount","MissingProp","HetCount","HetProp","MajorityBarcode","MajorityMeanFreq","Barcode"};
		
		OutputTextStore outStore = new UncompressedOutputTextStore(outRootFolder, "AllSamples.barcodes.tab");
		TableOutput out = new TableOutput (outStore, colNames, 1024 * 1024);
		for (int i = 0; i < samples.length; i++) {
			if (allGenos[i] == null) {
				continue;
			}
			out.newRow();
			out.appendValue(samples[i].getName());
			Genotype[] genos = allGenos[i];
			int bcodeLen = genos.length;
			
			int missingCount = 0;
			int hetCount = 0;
			double hetMajorFreqCum = 0.0;
			StringBuffer sbMaj = new StringBuffer(bcodeLen);
			StringBuffer sb = new StringBuffer(bcodeLen);
			for (int gIdx = 0; gIdx < bcodeLen; gIdx++) {
				sbMaj.append(genos[gIdx].topAllele);
				sb.append(genos[gIdx].genotype);

				if (genos[gIdx].call == 0) {
					missingCount++;
				} else if (genos[gIdx].call == 3) {
					double nraf = genos[gIdx].nraf;
					hetMajorFreqCum += (nraf > 0.5) ? nraf : (1.0 - nraf);
					hetCount++;
				}
			}
			double missingProp = ((double)missingCount) / ((double)bcodeLen);
			double hetProp = ((double)hetCount) / ((double)(bcodeLen-missingCount));

			out.appendValue(missingCount);
			out.appendValue(missingProp);
			out.appendValue(hetCount);
			out.appendValue(hetProp);
			out.appendValue(sbMaj.toString());
			if (hetCount > 0) {
				out.appendValue(hetMajorFreqCum / ((double)hetCount));
			} else {
				out.appendBlankValue();
			}
			out.appendValue(sb.toString());
		}
		out.close();
	}
	
	private void outputAlleleCounts (Sample[] samples, Genotype[][] allGenos) throws AnalysisException {
		String[] colNames = new String[snps.length + 1];
		colNames[0] = "Sample";
		for (int sIdx = 0; sIdx < snps.length; sIdx++) {
			colNames[sIdx+1] = snps[sIdx].getName();
		}
		OutputTextStore outStore = new UncompressedOutputTextStore(outRootFolder, "AllSamples.alleleCounts.tab");
		TableOutput out = new TableOutput (outStore, colNames, 1024 * 1024);
		for (int i = 0; i < samples.length; i++) {
			if (allGenos[i] == null) {
				continue;
			}
			out.newRow();
			out.appendValue(samples[i].getName());
			for (int sIdx = 0; sIdx < snps.length; sIdx++) {
				out.appendValue(allGenos[i][sIdx].alleleCountStr);
			}
		}
		out.close();
	}

	private void outputNrafs (Sample[] samples, Genotype[][] allGenos) throws AnalysisException {
		String[] colNames = new String[snps.length + 1];
		colNames[0] = "Sample";
		for (int sIdx = 0; sIdx < snps.length; sIdx++) {
			colNames[sIdx+1] = snps[sIdx].getName();
		}
		OutputTextStore outStore = new UncompressedOutputTextStore(outRootFolder, "AllSamples.nrafs.tab");
		TableOutput out = new TableOutput (outStore, colNames, 1024 * 1024);
		for (int i = 0; i < samples.length; i++) {
			if (allGenos[i] == null) {
				continue;
			}
			out.newRow();
			out.appendValue(samples[i].getName());
			for (int sIdx = 0; sIdx < snps.length; sIdx++) {
				out.appendValue(allGenos[i][sIdx].nraf);
			}
		}
		out.close();
	}

	public Genotype[] loadSampleResults (File resultFile) throws AnalysisException {
		Genotype[] genos = new Genotype[snps.length];
		String[] colNames = new String[] {"Chr","Pos","Genotype","TopAllele","Call","Nraf","Ref","Nref","Counts"};
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(resultFile));
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		int idx = 0;
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			
			//String chrName = values[0];
			//int pos = Integer.parseInt(values[1]); // TODO test position
			int call = Integer.parseInt(values[4]);
			char topAllele = values[3].charAt(0);
			int refCount = (call == 0) ? 0 : Integer.parseInt(values[6]);
			int nrefCount = (call == 0) ? 0 : Integer.parseInt(values[7]);
			String alleleCountStr = values[8];
			
			genos[idx++] = new Genotype(call, topAllele, refCount, nrefCount, alleleCountStr);
		}
		cfr.close();
		return genos;
	}

	/* ==========================================================
	 * Ancillary files
	 * ==========================================================
	 */
	public GenomePosition[] loadPositions (File snpListFile) throws AnalysisException {
		if (!snpListFile.canRead()) {
			throw new AnalysisException("Could not read SNP list file: " +snpListFile.getAbsolutePath());
		}
		ArrayList<GenomePosition> posList = new ArrayList<GenomePosition>();
		String[] colNames = new String[] {"Chr","Pos"};
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(snpListFile));
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			String chrName = values[0];
			int pos = Integer.parseInt(values[1]);
			posList.add(new GenomePosition(chrName,pos));
		}
		cfr.close();
		GenomePosition[] result = new GenomePosition[posList.size()];
		return posList.toArray(result);
	}

	/* ==========================================================
	 * Invocation
	 * ==========================================================
	 */
	public static class SingleSample {
		public static void main(String[] args) {
			if (args.length < 3) {
				System.err.println("Usage: org.cggh.bam.barcode.BarcodeFromVcfAnalysis$SingleSample <sampleName> <outputFolder> <snpListFile>");
				return;
			}
			String sampleName = args[0];	            log.info("Sample: "+sampleName);
			File outRootFolder = new File(args[1]);		log.info("Output Root Folder: "+outRootFolder.getAbsolutePath());
			File snpListFile = new File(args[2]);		log.info("SNP List File: "+snpListFile.getAbsolutePath());
			try {
				BarcodeFromVcfAnalysis task = new BarcodeFromVcfAnalysis(outRootFolder, snpListFile);
				Sample sample = new Sample (sampleName, null, null);
				task.analyzeSample(sample);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
	
	public static class MultiSample {
		public static void main(String[] args) {
			if (args.length < 3) {
				System.err.println("Usage: org.cggh.bam.barcode.BarcodeFromVcfAnalysis$MultiSample <sampleListFile> <outputFolder> <snpListFile>");
				return;
			}
			File sampleListFile = new File(args[0]);	log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File outRootFolder = new File(args[1]);		log.info("Output Root Folder: "+outRootFolder.getAbsolutePath());
			File snpListFile = new File(args[2]);		log.info("SNP List File: "+snpListFile.getAbsolutePath());
			
			int maxThreads = Integer.parseInt(System.getProperty("maxThreads", "0"));
			
			try {
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				BarcodeFromVcfAnalysis task = new BarcodeFromVcfAnalysis(outRootFolder, snpListFile);
				multi.execute((SampleAnalysis) task);
				task.analyzeAllSampleResults(multi.getSamples());
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}

	
	public static class MergeResults {
		public static void main(String[] args) {
			if (args.length < 3) {
				System.err.println("Usage: org.cggh.bam.barcode.BarcodeFromVcfAnalysis$MergeResults <sampleListFile> <outputFolder> <snpListFile>");
				return;
			}
			File sampleListFile = new File(args[0]);	log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File outRootFolder = new File(args[1]);		log.info("Output Root Folder: "+outRootFolder.getAbsolutePath());
			File snpListFile = new File(args[2]);		log.info("SNP List File: "+snpListFile.getAbsolutePath());
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				BarcodeFromVcfAnalysis task = new BarcodeFromVcfAnalysis(outRootFolder, snpListFile);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}

}
