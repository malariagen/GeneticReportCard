package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.GenomeRegion;
import org.cggh.common.textStore.*;
import org.cggh.common.threading.ParallelExecutableManager;
import org.cggh.common.util.FileUtilities;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;


public class ReadsAlignmentAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private Sample[]      samples;
	private GenotypingConfig config;
	private Locus[]       loci;
	private File          outLociFolder;
	private boolean       isSingleSample;
	
	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	
	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public ReadsAlignmentAnalysis (File configFile, String sampleId, File sampleBamFile, File refFastaFile, File outRootFolder) throws AnalysisException  {
		this(configFile, new Sample[] { new Sample (sampleId, sampleBamFile) }, refFastaFile, outRootFolder);
		isSingleSample = true;
	}
	
	/* ==========================================================
	 * Invocation: sample/bamFile list file
	 * ==========================================================
	 */
	public ReadsAlignmentAnalysis (File configFile, File sampleListFile, File refFastaFile, File outRootFolder) throws AnalysisException  {
		this(configFile, readSamples(sampleListFile), refFastaFile, outRootFolder);
		isSingleSample = false;
	}

	/* ==========================================================
	 * Initialization
	 * ==========================================================
	 */
	private ReadsAlignmentAnalysis (File configFile, Sample[] samples, File refFastaFile, File outRootFolder) throws AnalysisException  {

		// Get the list of samples
		this.samples = samples;

		// Load up the reference genome sequences
		ReferenceGenome.initialize(refFastaFile);
		
		// Perse configuration file
		config = new GenotypingConfig (configFile);
		loci = config.getLoci();
		
		// Create the root output directory
		if (!outRootFolder.exists()) {
			outRootFolder.mkdirs();
		}
		this.outLociFolder = FileUtilities.checkFolder(outRootFolder, "loci", true);		
	}

	
	public void execute () throws AnalysisException, IOException  {
		if (isSingleSample) {
			executeSingleSample();
		} else {
			executeMultiSample();
		}
	}
	
	public void executeSingleSample () throws AnalysisException, IOException  {
		// Create one task per sample
		SampleAnalysisTask[] analysisTasks = new SampleAnalysisTask[1];
		analysisTasks[0] = new SampleAnalysisTask (0);
		analysisTasks[0].run();
	}
	
	public void executeMultiSample () throws AnalysisException, IOException  {

		// Analyze the samples one by one at all the loci		
		ParallelExecutableManager pem;
		int threadCount = config.getMaxThreads();
		if (threadCount > 0) { 
			pem = new ParallelExecutableManager(threadCount);
		} else {
			pem = new ParallelExecutableManager();
		}
		
		// Create one task per sample
		SampleAnalysisTask[] analysisTasks = new SampleAnalysisTask[samples.length];
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			analysisTasks[sIdx] = new SampleAnalysisTask (sIdx);
			pem.addTask(analysisTasks[sIdx]);
		}

		// Analyze all the samples at all the loci
		pem.startExecution();
		pem.waitForCompletion();
	}
	
	
	private class SampleAnalysisTask implements Runnable {
		int sampleIdx;

		public SampleAnalysisTask (int sampleIdx) {
			this.sampleIdx = sampleIdx;
		}
		
		@Override
		public void run() {
			log.info("Starting " + samples[sampleIdx].getName() + " (" + (sampleIdx+1) + " of " + samples.length + ")");
			try {
				analyzeSample (sampleIdx);
			} catch (Exception e) {
				log.info("Aborting " + samples[sampleIdx].getName());
				log.error("Error processing BAM file for sample "+ samples[sampleIdx].getName() + ": "+e);
				e.printStackTrace();
			}
			log.info("Completed " + samples[sampleIdx].getName());
		}
	}

	
	@SuppressWarnings("unchecked")
	public void analyzeSample (int sampleIdx) throws AnalysisException, IOException  {
		
		// **************************************
		// PART 1 - Reads analysis and Alignment
		// **************************************
		
		// Read the reads from the SAM file
		Sample sample = samples[sampleIdx];
		ArrayList<MappedRead>[] mappedReadLists = new ArrayList[loci.length];
		
		for (int i = 0; i < loci.length; i++) {
			Locus locus = loci[i];
			mappedReadLists[i] = new ArrayList<MappedRead>();
			getMappedLocusReads (sample, locus, mappedReadLists[i]);
		}
		
		// Then search unmapped reads (this does not necessarily find them all, but will do)
		getUnmappedLocusReads (sample, loci, mappedReadLists);
		
		// Perform a reads alignment and count the alleles in the specified region
		for (int i = 0; i < loci.length; i++) {
			Locus locusConfig = loci[i];
			//log.info("Locus: "+locusResult.config.name+" - Targets: "+locusConfig.targets.length);
				
			// Get the list of reads for this genotype(locus region)
			Iterator<MappedRead> it = mappedReadLists[i].iterator();
			while(it.hasNext()) {
				if (!it.next().coversAtLeastOneTargets()) {
					it.remove();
				}				
			}
			MappedRead[] mappedReads = mappedReadLists[i].toArray(new MappedRead[mappedReadLists[i].size()]);
			//log.info("Got reads: "+mappedReads.length);
			
			// Construct a read alignment
			ReadsAlignment ra = new ReadsAlignment(mappedReads, locusConfig, sample);
			//log.info("Misaligned: "+ra.misalignCount);
			
			// Write out to file the reads alignment
			ra.outputAlignment(outLociFolder);
			//log.info("Written alignments");
		}
	}
	
	private void getMappedLocusReads (Sample sample, Locus locus, ArrayList<MappedRead> mappedReadsList) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		GenomeRegion readSearchInterval =  locus.getReadSearchInterval();
		Target[] targets = locus.getTargets();
		SAMRecordIterator it = samReader.query(readSearchInterval.getChromosome(), readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true);
		while (it.hasNext()) {
			SAMRecord record = it.next();
			boolean matched = false;
			
			// If the read is ungapped, we can use the BAM start/end coords
			if ((record.getCigarLength() == 1) && (record.getCigarString().endsWith("M"))) {
				int startPos = record.getAlignmentStart();
				int endPos = record.getAlignmentEnd();
				
				// First test: does it cover at least one genotyping target?
				for (int i = 0; i < targets.length; i++) {
					GenomeRegion targetRegion = targets[i].getTargetRegion();
					if ((startPos <= targetRegion.getStartPos()) && (endPos >= targetRegion.getStopPos())) {
						MappedRead sr = new MappedRead(record, locus, startPos, MappedRead.MAPPED);
						mappedReadsList.add(sr);
						matched = true;
						break;
					}
				}
				if (matched) {
					continue;
				}
			}
			
			// If not, take the ungapped read, andd treat it as if unmapped, try to find an anchor
			matched = matchUnmappedReadAtLocus (record, locus, mappedReadsList, MappedRead.REMAPPED);
		}
		it.close();
	}


	private void getUnmappedLocusReads (Sample sample, Locus[] loci, ArrayList<MappedRead>[] mappedReadLists) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		SAMRecordIterator it = samReader.queryUnmapped();
		while (it.hasNext()) {
			SAMRecord record = it.next();
			boolean matched = matchUnmappedRead (record, loci, mappedReadLists);
					 
			if (!matched) {
				SAMRecordUtil.reverseComplement(record);
				matched = matchUnmappedRead (record, loci, mappedReadLists);
			}
		}
		it.close();
	}

	private boolean matchUnmappedRead (SAMRecord record, Locus[] loci, ArrayList<MappedRead>[] mappedReadLists) throws AnalysisException {
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus locus = loci[lIdx];
			ArrayList<MappedRead> mappedReadList = mappedReadLists[lIdx];
			if (matchUnmappedReadAtLocus (record, locus, mappedReadList, MappedRead.UNMAPPED)) {
				return true;
			}
		}
		return false;
	}

	private boolean matchUnmappedReadAtLocus (SAMRecord record, Locus locus, ArrayList<MappedRead> mappedReadList, int mappingStatus) throws AnalysisException {

		// Does the read contain an anchor?
		String readSequence = record.getReadString();
		Anchor[] anchors = locus.getAnchors();
		for (int aIdx = 0; aIdx < anchors.length; aIdx++) {
			Matcher m = anchors[aIdx].getRegex().matcher(readSequence);
		    if (m.find()) {
		    	int anchorPos = m.start();
				MappedRead sr = new MappedRead(record, locus, anchors[aIdx], anchorPos, mappingStatus);
				mappedReadList.add(sr);
				return true;
		    }									
		}					
		return false;
	}

	
	/* *************************************************************************
	 * Initialization routines
	 * *************************************************************************
	 */
	private static Sample[] readSamples(File sampleListFile) throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));
		String[] colNames = new String[]{"Sample","BamFile"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample(values[0], new File(values[1]));
			if (!s.getBamFile().canRead()) {
				throw new AnalysisException("Error getting samples list: BAM file " + values[1] + " cannot be read.");
			}
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}




	/* ==========================================================
	 * Execution
	 * ==========================================================
	 */
	public static void main(String[] args) {
		if (args.length < 4) {
			showUsage();
			return;
		}
		
		ReadsAlignmentAnalysis task = null;
		try {
			String mode = args[0];
			if ("singleSample".equals(mode)) {
				String sampleName = args[1];
				File sampleBamFile = new File(args[2]);
				File configFile = new File(args[3]);
				File refFastaFile = new File(args[4]);
				File rootFolder = new File(args[5]);
				task = new ReadsAlignmentAnalysis(configFile, sampleName, sampleBamFile, refFastaFile, rootFolder);
				
			} else if ("sampleList".equals(mode)) {
				File sampleListFile = new File(args[1]);
				File configFile = new File(args[2]);
				File refFastaFile = new File(args[3]);
				File rootFolder = new File(args[4]);
				task = new ReadsAlignmentAnalysis(configFile, sampleListFile, refFastaFile, rootFolder);

			} else {
				showUsage();
				return;
			}
			
			task.execute();
			
		} catch (Exception e) {
			log.error("Error executing task: " + e);
			return;
		}
		log.info("Exiting");
	}
	
	public static void showUsage() {
		System.out.println("Usage: org.cggh.bam.genotyping.reads.ReadsAlignmentAnalysis [singleSample|sampleList] [<sampleName> <bamFile>|<sampleListFile>] <configFile> <refFasta> <rootFolder>");
	}
	

}
