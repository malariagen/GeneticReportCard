package org.cggh.bam.readCounts;

import org.cggh.bam.*;
import org.cggh.bam.genotyping.*;
import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import htsjdk.samtools.*;
import java.io.*;
import java.text.*;
import java.util.*;
import org.apache.commons.logging.*;


public class ReadCountAnalysis extends SampleAnalysis {
	
	private static Log log = LogFactory.getLog((String)ClassUtilities.getCurrentClassName());
	
	public static final int MIN_PHRED_SCORE = 20;
	private static final SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	private GenotypableSnp[] genoPositions;
	private ReadCountConfig config;
	private String prefix = null;
	private File warningFile = null;
	
	public ReadCountAnalysis(File configFile, File snpListFile, File outRootFolder, String prefix) throws AnalysisException {
		super(configFile, outRootFolder);
		this.config = new ReadCountConfig (configFile);
		this.prefix = prefix;
		genoPositions = readGenotypableSnp(snpListFile);
	}

	public void analyzeSample(Sample sample) throws AnalysisException {
		try {
			analyzeSampleReads(sample);
		} catch (AnalysisException e) {
	    	recordWarning (sample, null, "Error analyzing sample: "+e, true);
	    	throw e;
		}
	}
	
	public void analyzeSampleReads(Sample sample) throws AnalysisException {
			
	    // See GentypingConfig for parameters. Right now, it's min 5 reads for a call, min 5% total reads to call an allele
	    NucleotideGenotyper bg = new NucleotideGenotyper(config);

	    log.info("Starting " + sample.getName());  
		int[] refCounts = new int[genoPositions.length];
		int[] nrefCounts = new int[genoPositions.length];
		double[] genoFreq = new double[genoPositions.length];  Arrays.fill(genoFreq, Double.NaN);
		int[] genoNum = new int[genoPositions.length];
		int[] genoNumMulti = new int[genoPositions.length];
		File bamFile = sample.getBamFile();
		if ((!bamFile.exists()) || (!bamFile.canRead())) {
	    	recordWarning (sample, null, "Cannot read bam file "+bamFile.getAbsolutePath(), true);
	    	return;
		}
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		for (int snpIdx = 0; snpIdx <  genoPositions.length; snpIdx++) {
			GenotypableSnp genoPos = genoPositions[snpIdx];
			
			// Get the read counts by inspecting the BAM
		    NtAlleleCounter ntCounts = getReadCounts (samReader, genoPos);
		    refCounts[snpIdx] = ntCounts.getCount(genoPos.ref);
		    nrefCounts[snpIdx] = ntCounts.getCount(genoPos.nonref);
		    if (!bg.hasSufficientReads(ntCounts)) {
		    	continue;
		    }
		    // We can genotype. First check if there is an unexpected majority allele
		    char maj = bg.getMajorityAllele(ntCounts);
		    if ((maj != genoPos.ref) && (maj != genoPos.nonref) && (maj != '-')) {
		    	recordWarning (sample, genoPos, 
		    			"Majority allele "+ maj+ " ("+ ntCounts.getCount(maj)+ " out of "+ ntCounts.getCumulativeCount()+ " reads) "
		    			+ "is not ref/nref ["+ genoPos.ref+ "/"+ genoPos.nonref+ "]", false);
		    	continue;
		    }
		    double freqCall = bg.getFrequencyCall(ntCounts, genoPos.ref, genoPos.nonref);
		    genoFreq[snpIdx] = freqCall;
		    genoNum[snpIdx] = bg.frequencyToMajorityAlleleNum(freqCall);
		    genoNumMulti[snpIdx] = bg.frequencyToAlleleNum(freqCall);
		}
		
		outputReadCounts (sample, refCounts, nrefCounts, genoFreq, genoNum, genoNumMulti);
		log.info("Completed " + sample.getName());
	}
	
	private NtAlleleCounter getReadCounts (SamReader samReader, GenotypableSnp genoPos) {
		int pos = genoPos.getPos();
		NtAlleleCounter counter = new NtAlleleCounter();
		SAMRecordIterator rIt = samReader.query(genoPos.getChromosome(), pos, pos, false);
		while (rIt.hasNext()) {
			final SAMRecord rec = (SAMRecord) rIt.next();
			Iterator<AlignmentBlock> aIt = rec.getAlignmentBlocks().iterator();
			while (aIt.hasNext()) {
				AlignmentBlock ab = aIt.next();
				int abStart = ab.getReferenceStart();
				if (pos < abStart) {
					break;
				}
				int abEnd = abStart + ab.getLength() - 1;
				if (pos <= abEnd) {
					int offset = pos - abStart;
					int rOffset = offset + ab.getReadStart() - 1;
					
					String readQ = rec.getBaseQualityString();
					byte baseQ = (byte)(readQ.charAt(rOffset) - 33);
					if (baseQ < MIN_PHRED_SCORE) {
						break;
					}
					String read = rec.getReadString();
					counter.increment(read.charAt(rOffset));
				}
			}
		}
		rIt.close();
		return counter;
	}
	
	private void outputReadCounts (Sample sample, int[] refCounts, int[] nrefCounts, double[] genoFreq, int[] genoNum, int[] genoNumMulti) throws AnalysisException {
		File sampleFile = getSampleResultsFile (sample);
		TableOutput out = new TableOutput(sampleFile.getParentFile(), sampleFile.getName(), READCOUNT_HEADERS, 64*1024);
		for (int snpIdx = 0; snpIdx <  genoPositions.length; snpIdx++) {
			out.newRow();
			out.appendValue(genoPositions[snpIdx].getChromosome());
			out.appendValue(genoPositions[snpIdx].getPos());
			out.appendValue(refCounts[snpIdx]);
			out.appendValue(nrefCounts[snpIdx]);
			out.appendValue(genoFreq[snpIdx]);
			out.appendValue(genoNumMulti[snpIdx]);
			out.appendValue(genoNum[snpIdx]);
		}
		out.close();
	}

	private File getSampleResultsFile (Sample sample) {
		File sampleFolder = getSampleSubfolder(outRootFolder, sample, true);
		File sampleFile = new File(sampleFolder, sample.getName()+'.'+prefix+".tab");
		return sampleFile;
	}

	public GenotypableSnp[] readGenotypableSnp (File snpListFile) throws AnalysisException {
		if (!snpListFile.canRead()) {
			throw new AnalysisException("Could not read SNP list file: " +snpListFile.getAbsolutePath());
		}
		ArrayList<GenotypableSnp> posList = new ArrayList<GenotypableSnp>();
		String[] colNames = new String[] {"Chr","Pos","Ref","Nonref"};
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(snpListFile));
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			String chrName = values[0];
			int pos = Integer.parseInt(values[1]);
			char ref = values[2].trim().charAt(0);
			char nref = values[3].trim().charAt(0);
			posList.add(new GenotypableSnp(chrName,pos,ref,nref));
		}
		cfr.close();
		Collections.sort(posList);
		GenotypableSnp[] result = new GenotypableSnp[posList.size()];
		return posList.toArray(result);
	}
	
	private void recordWarning (Sample sample, GenotypableSnp genoPos, String warning, boolean isFatal) throws AnalysisException {
    	if (warningFile == null) {
        	String warningFilename = prefix+".warnings.tab";
        	warningFile = new File (outRootFolder, warningFilename);
    	}
    	String posStr = genoPos == null ? "<none>" : genoPos.getName();
    	String msg = '\n'+df.format(new Date())+'\t'+sample.getName()+'\t'+posStr+'\t'+warning;
		try {
    	    if (!warningFile.exists()) {
				FileUtilities.writeFileContent("Timestamp\tSample\tPos\tMessage", warningFile);
        	}
    	    FileUtilities.appendFileContent(msg, warningFile);
		} catch (IOException e) {
			throw new AnalysisException("Error writing warning to "+ warningFile.getAbsolutePath() + ": "+e);
		}
		String displayMsg = "Sample " + sample.getName() + ((genoPos==null)?"":" - Pos "+posStr) + " - " + msg;
		log.warn(displayMsg);
		if (isFatal) {
			log.error("Aborting " + sample.getName());
		}
	}

	/* ==========================================================
	 * Result Merging
	 * ==========================================================
	 */
	private static final String[] READCOUNT_HEADERS = new String[] { "Chr", "Pos", "Ref", "Nonref", "GenotypeFreq", "GenotypeNumMulti", "GenotypeNum" };
	
	public void mergeAllSampleResults(Sample[] samples) throws AnalysisException, IOException {

		ArrayList<String> sampleNameList = new ArrayList<String>();
		ArrayList<File> sampleFileList = new ArrayList<File>();
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFile = getSampleResultsFile (sample);
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
				continue;
			}
			sampleNameList.add(sample.getName());
			sampleFileList.add(sampleFile);
		}
		String[] sampleNames = sampleNameList.toArray(new String[sampleNameList.size()]);
		File[] sampleFiles = sampleFileList.toArray(new File[sampleFileList.size()]);
		
		short[][] refCounts   = new short[genoPositions.length][sampleNames.length];	
		short[][] nrefCounts  = new short[genoPositions.length][sampleNames.length];
		float[][] genoFreq    = new float[genoPositions.length][sampleNames.length];
		byte[][] genoNumMulti = new byte[genoPositions.length][sampleNames.length];
		byte[][] genoNum      = new byte[genoPositions.length][sampleNames.length];

		for (int sIdx = 0; sIdx < sampleFiles.length; sIdx++) {
			ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleFiles[sIdx]));
			ColumnFileReader.ColumnReader cr = cfr.getColumnReader(READCOUNT_HEADERS);
			int gpIdx = 0;
			while (cfr.nextRecord()) {
				String[] values = cr.getValues();
				String chrName = values[0];
				int pos = Integer.parseInt(values[1]);
				GenomePosition gp = genoPositions[gpIdx];
				if (!chrName.equals(gp.getChromosome()) || pos != gp.getPos()) {
					throw new AnalysisException ("Error aggregating data from "	+ sampleFiles[sIdx].getAbsolutePath()
							+ " at line " + cfr.getLineNumber()+ ": found position " + chrName + ":" + pos+ " instead of "+ gp.getName());
				}
				refCounts[gpIdx][sIdx]    = Short.parseShort(values[2]);
				nrefCounts[gpIdx][sIdx]   = Short.parseShort(values[3]);
				genoFreq[gpIdx][sIdx]     = (values[4].equals("-")) ? Float.NaN : Float.parseFloat(values[4]);
				genoNumMulti[gpIdx][sIdx] = Byte.parseByte(values[5]);
				genoNum[gpIdx][sIdx]      = Byte.parseByte(values[6]);
				gpIdx++;
			}
			cfr.close();
		}
		outputMergedValues (sampleNames, "ReadCounts.ref",   OUT_TYPE_SHORT, refCounts);
		outputMergedValues (sampleNames, "ReadCounts.nref",  OUT_TYPE_SHORT, nrefCounts);
		outputMergedValues (sampleNames, "GenotypeFreq",     OUT_TYPE_FLOAT, genoFreq);
		outputMergedValues (sampleNames, "GenotypeNumMulti", OUT_TYPE_BYTE,  genoNumMulti);
		outputMergedValues (sampleNames, "GenotypeNum",      OUT_TYPE_BYTE,  genoNum);
	}
	
	private static final int OUT_TYPE_SHORT = 1;
	private static final int OUT_TYPE_FLOAT = 2;
	private static final int OUT_TYPE_BYTE = 3;
	
	private void outputMergedValues (String[] sampleNames, String fileTypePart, int type, Object[] valueArray) throws AnalysisException {
		String[] colNames = TextUtilities.mergeStringLists(new String[] {"Chr", "Pos"}, sampleNames);	
		OutputTextStore ots = new GzippedOutputTextStore (outRootFolder, prefix+'.'+fileTypePart+".tab");
		TableOutput out = new TableOutput(ots, colNames, 64*1024);
		for (int snpIdx = 0; snpIdx <  genoPositions.length; snpIdx++) {
			out.newRow();
			out.appendValue(genoPositions[snpIdx].getChromosome());
			out.appendValue(genoPositions[snpIdx].getPos());
			switch (type) {
			case OUT_TYPE_SHORT:
				short[] counts = (short[])valueArray[snpIdx];
				for (int i = 0; i < sampleNames.length; i++) {
					out.appendValue(counts[i]);
				}
				break;
			case OUT_TYPE_FLOAT:
				float[] freqs = (float[])valueArray[snpIdx];
				for (int i = 0; i < sampleNames.length; i++) {
					out.appendValue(freqs[i]);
				}
				break;
			case OUT_TYPE_BYTE:
				byte[] genos = (byte[])valueArray[snpIdx];
				for (int i = 0; i < sampleNames.length; i++) {
					out.appendValue(genos[i]);
				}
				break;
			}
		}
		out.close();
	}


	/* ==========================================================
	 * Data Structures
	 * ==========================================================
	 */
	private class GenotypableSnp extends GenomePosition {
		char ref;
		char nonref;
		public GenotypableSnp(String chr, int pos, char ref, char nonref) {
			super(chr, pos);
			this.ref = ref;
			this.nonref = nonref;
		}
	}
	

	/* ==========================================================
	 * Single Sample Execution
	 * ==========================================================
	 */
	public static class SingleSample {
		public static void main(String[] args) {
			if (args.length < 7) {
				log.error("Usage: org.cggh.bam.readCounts.ReadCountAnalysis$SingleSample <configFile> <prefix> <batchId> <sampleName> <bamFile> <snpListFile> <outFolder>");
				return;
			}
			File configFile = new File(args[0]);		    log.info("ConfigFile: "+configFile.getAbsolutePath());
			String prefix = args[1].trim();		            log.info("Prefix: " + prefix);
			String batchId = args[2];						log.info("BatchId: "+batchId);
			String sampleId = args[3];						log.info("SampleId: " + sampleId);
			File sampleBamFile = new File(args[4]);			log.info("SampleBamFile: " + sampleBamFile.getAbsolutePath());
			File snpListFile = new File(args[5]);			log.info("SnpListFile: " + snpListFile.getAbsolutePath());
			File outRootFolder = new File(args[6]);	        log.info("OutRootFolder: " + outRootFolder.getAbsolutePath());
			try {
				Sample sample = new Sample(batchId, sampleId, sampleBamFile);
				ReadCountAnalysis task = new ReadCountAnalysis(configFile, snpListFile, outRootFolder, prefix);
				task.analyzeSample(sample);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				e.printStackTrace();
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
				log.error("Usage: org.cggh.bam.readCounts.ReadCountAnalysis$MultiSample <configFile> <prefix> <batchId> <sampleListFile> <snpListFile> <rootFolder>");
				return;
			}
			File configFile = new File(args[0]);		    log.info("ConfigFile: "+configFile.getAbsolutePath());
			String prefix = args[1].trim();		            log.info("Prefix: " + prefix);
			File sampleListFile = new File(args[2]);		log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File snpListFile = new File(args[3]);			log.info("SnpListFile: " + snpListFile.getAbsolutePath());
			File outRootFolder = new File(args[4]);			log.info("OutRootFolder: " + outRootFolder.getAbsolutePath());
			int maxThreads = Integer.parseInt(System.getProperty("maxThreads", "0"));

			try {
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				ReadCountAnalysis task = new ReadCountAnalysis(configFile, snpListFile, outRootFolder, prefix);
				multi.execute((SampleAnalysis) task);
				task.mergeAllSampleResults(multi.getSamples());
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				e.printStackTrace();
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
				log.error("Usage: org.cggh.bam.readCounts.ReadCountAnalysis$MergeResults <configFile> <prefix> <sampleListFile> <snpListFile> <outRootFolder>");
				return;
			}
			File configFile = new File(args[0]);		    log.info("ConfigFile: "+configFile.getAbsolutePath());
			String prefix = args[1].trim();		            log.info("Prefix: " + prefix);
			File sampleListFile = new File(args[2]);		log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
			File snpListFile = new File(args[3]);			log.info("SnpListFile: " + snpListFile.getAbsolutePath());
			File outRootFolder = new File(args[4]);			log.info("OutRootFolder: " + outRootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				ReadCountAnalysis task = new ReadCountAnalysis(configFile, snpListFile, outRootFolder, prefix);
				task.mergeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				e.printStackTrace();
				return;
			}
			log.info("Exiting");
		}
	}
}
