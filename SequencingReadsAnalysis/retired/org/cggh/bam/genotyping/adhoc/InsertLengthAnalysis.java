package org.cggh.bam.genotyping.adhoc;

import org.cggh.bam.genotyping.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.GenomeRegion;
import org.cggh.common.textStore.*;
import org.cggh.common.threading.ParallelExecutableManager;
import org.cggh.common.util.Statistics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class InsertLengthAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private Sample[]      samples;
	private GenotypingConfig config;
	private Locus[] locusConfigs;
	private File          outRootFolder;
	
	private TableOutput meansOut = null;
	private TableOutput mediansOut = null;
	private TableOutput stdevsOut = null;
	private TableOutput countsOut = null;
	private TableOutput largeInsertOut = null;
	
	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	

	/* ==========================================================
	 * Samples and SNP lists initialization
	 * ==========================================================
	 */
	public InsertLengthAnalysis (File configFile, File sampleListFile, File refFastaFile, File rootFolder) throws AnalysisException  {

		// Load up the reference genome sequences
		ReferenceGenome.initialize(refFastaFile);
		
		// Perse configuration file
		config = new GenotypingConfig (configFile);
		locusConfigs = config.getLoci();
		
		// Get the list of samples
		samples = readSamples(sampleListFile);

		// Create the root output directory
		this.outRootFolder = rootFolder;
		if (!rootFolder.exists()) {
			rootFolder.mkdirs();
		}
	}

	
	public void execute () throws AnalysisException, IOException  {
	
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

		initializeOutputs ();

		// Analyze all the samples at all the loci
		pem.startExecution();
		pem.waitForCompletion();

		meansOut.close();
		mediansOut.close();
		stdevsOut.close();
		countsOut.close();
		largeInsertOut.close();

		return;
	}
	
	
	private class SampleAnalysisTask implements Runnable {
		int                 sampleIdx;
		SampleLocusResult[] results;

		public SampleAnalysisTask (int sampleIdx) {
			this.sampleIdx = sampleIdx;
			this.results = new SampleLocusResult[locusConfigs.length];
		}
		
		@Override
		public void run() {
			Sample sample = samples[sampleIdx];
			log.info("Starting " + sample.getName() + " (" + (sampleIdx+1) + " of " + samples.length + ")");
			try {
				// Read the reads from the SAM file
				for (int i = 0; i < locusConfigs.length; i++) {
					results[i] = processMappedLocusReads (sample, locusConfigs[i]);
				}
				outputSampleResults (sample, results);

			} catch (Exception e) {
				log.info("Aborting " + samples[sampleIdx].getName());
				log.error("Error processing BAM file for sample "+ samples[sampleIdx].getName() + ": "+e);
				e.printStackTrace();
			}
			log.info("Completed " + samples[sampleIdx].getName());
		}
	}

	public class SampleLocusResult {
		int count;
		double insertLenMean;
		double insertLenMedian;
		double insertLenSd;

		public SampleLocusResult(int[] lengths) {
			count = lengths.length;
			Statistics stat = new Statistics(lengths);
			insertLenMean = stat.getMean();
			insertLenMedian = stat.getMedian();
			insertLenSd = stat.getStdev();
		}
	};
	
	
	private SampleLocusResult processMappedLocusReads (Sample sample, Locus locusConfig) throws AnalysisException {
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		GenomeRegion readSearchInterval =  locusConfig.getReadSearchInterval();
		
		HashMap<String,Integer> lengthMap = new HashMap<String,Integer>();
		SAMRecordIterator it = samReader.query(readSearchInterval.getChromosome(), readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true);
		while (it.hasNext()) {
			SAMRecord record = it.next();
			String id = record.getReadName();
			if (lengthMap.containsKey(id)) {
				continue; // Mate already done
			}
			int insertLength = record.getInferredInsertSize();
			if (insertLength == 0) {
				continue;
			}
			insertLength = (insertLength < 0) ? -insertLength : insertLength;
			if (insertLength > 1000) {
				outputLargeInsertInfo (locusConfig, sample, record);
				continue;	// Forget mismappings or major rearrangements
			}
			lengthMap.put(id, insertLength);
		}
		it.close();
		
		// Process the insert length data
		int cnt = lengthMap.size();
		int[] lengths = new int[cnt];
		List<Integer> lengthList = new ArrayList<Integer>(lengthMap.values());
		for (int i = 0; i < cnt; i++) {
			lengths[i] = lengthList.get(i);
		}
		SampleLocusResult result = new SampleLocusResult(lengths);
		return result;
	}


	private void initializeOutputs () throws AnalysisException {
		// Get the file headers
		ArrayList<String> colNameList = new ArrayList<String>();
		colNameList.add("Sample");
		for (int lIdx = 0; lIdx < locusConfigs.length; lIdx++) {
			colNameList.add(locusConfigs[lIdx].getName());
		}
		String[] colNames = colNameList.toArray(new String[colNameList.size()]);
		meansOut = new TableOutput (outRootFolder, "InsertLengthMeans.tab", colNames, 128*1024);
		mediansOut = new TableOutput (outRootFolder, "InsertLengthMedians.tab", colNames, 128*1024);
		stdevsOut = new TableOutput (outRootFolder, "InsertLengthStdevs.tab", colNames, 128*1024);
		countsOut = new TableOutput (outRootFolder, "InsertLengthCounts.tab", colNames, 128*1024);
		
		colNames = new String[] {"Sample","Locus","InsertSize","ReadId"};
		largeInsertOut = new TableOutput (outRootFolder, "LargeInsertReads.tab", colNames, 1024*1024);			
	}


	private synchronized void outputSampleResults (Sample sample, SampleLocusResult[] sampleResults) throws AnalysisException {
		meansOut.newRow();		meansOut.appendValue(sample.getName());
		mediansOut.newRow();	mediansOut.appendValue(sample.getName());
		stdevsOut.newRow();		stdevsOut.appendValue(sample.getName());
		countsOut.newRow();		countsOut.appendValue(sample.getName());
		for (int lIdx = 0; lIdx < sampleResults.length; lIdx++) {
			SampleLocusResult sResult = sampleResults[lIdx];
			meansOut.appendValue(sResult.insertLenMean);
			mediansOut.appendValue(sResult.insertLenMedian);
			stdevsOut.appendValue(sResult.insertLenSd);
			countsOut.appendValue(sResult.count);
		}
		meansOut.commitIfHasContent();
		mediansOut.commitIfHasContent();
		stdevsOut.commitIfHasContent();
		countsOut.commitIfHasContent();
	}
	
	private synchronized void outputLargeInsertInfo (Locus locusConfig, Sample sample, SAMRecord record) throws AnalysisException {
		//log.info("Found long insert length at locus "+ locusConfig.getName()+" for sample "+ sample.getName()+ ": "+ record.getInferredInsertSize()+ " (read ID: "+ record.getReadName()+ ")");
		largeInsertOut.newRow();
		largeInsertOut.appendValue(sample.getName());
		largeInsertOut.appendValue(locusConfig.getName());
		largeInsertOut.appendValue(record.getInferredInsertSize());
		largeInsertOut.appendValue(record.getReadName());
	}
	

	/* *************************************************************************
	 * Initialization routines
	 * *************************************************************************
	 */
	private Sample[] readSamples(File sampleListFile) throws AnalysisException {
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
			System.out.println("Usage: org.cggh.bam.genotyping.adhoc.InsertLengthAnalysis <configFile> <sampleListFile> <refFasta> <rootFolder>");
			return;
		}
		File configFile = new File(args[0]);
		File sampleListFile = new File(args[1]);
		File refFastaFile = new File(args[2]);
		File rootFolder = new File(args[3]);
		
		try {
			InsertLengthAnalysis task = new InsertLengthAnalysis(configFile, sampleListFile, refFastaFile, rootFolder);
			task.execute();
		} catch (Exception e) {
			log.error("Error: " + e);
			return;
		}
		log.info("Exiting");
	}

}
