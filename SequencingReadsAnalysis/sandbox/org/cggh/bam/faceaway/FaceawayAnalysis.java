package org.cggh.bam.faceaway;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

//import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class FaceawayAnalysis extends SampleAnalysis  {
	
	//private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	private SamReaderFactory  samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	public FaceawayAnalysis (File outRootFolder) throws AnalysisException  {
		super (outRootFolder);
	}

	/* **********************************************************************
	 * Single sample processing
	 * **********************************************************************
	 */
	public void analyzeSample (Sample sample) throws AnalysisException  {
		log.info("Starting " + sample.getName());
		try {
			// First reader: for scanning for faceaways
			SamReader reader1 = samReaderFactory.open(sample.getBamFile());
			SAMRecordIterator it = reader1.iterator();
			
			// Second reader: for finding mates
			SamReader reader2 = samReaderFactory.open(sample.getBamFile());
			
			//
			ArrayList<Faceaway> fList = new ArrayList<Faceaway>();
			while (it.hasNext()) {
				SAMRecord record = it.next();
				// This will ignore all secondary and supplementary alignments, reads not passing vendor filters, as well as any duplicates
				if (record.getFlags() >= 256) {
					continue;
				}
				
				// Check for faceaway
				if (record.getReadNegativeStrandFlag() && !(record.getMateNegativeStrandFlag())
					&& !(record.getNotPrimaryAlignmentFlag())
					&& record.getReferenceName().equals(record.getMateReferenceName())
					&& (record.getMateAlignmentStart() > record.getAlignmentEnd())
					&& (record.getMappingQuality() > 0)) {
					
					// We have the negative strand, find the mate and check for quality
					SAMRecord mate = reader2.queryMate(record);
					if ((mate == null) 
						|| (mate.getNotPrimaryAlignmentFlag())
						|| (mate.getMappingQuality() <= 0)) {
						continue;
					}
					
					Faceaway f = new Faceaway(record, mate);
					fList.add(f);
				}
			}
			it.close();
			reader1.close();
			reader2.close();
			
			Faceaway[] fResults = fList.toArray(new Faceaway[fList.size()]);
			log.info("Found faceaways: " + fResults.length);
			SampleResult sr = new SampleResult(sample, fResults);

			// Write out the results
			File outFolder = getSampleSubfolder (outRootFolder, sample, true);
			sr.writetoFile(outFolder);
			
		} catch (Exception e) {
			log.info("Aborting " + sample.getName());
			log.error("Error processing BAM file for sample "+ sample.getName() + ": "+e);
			e.printStackTrace();
		}
		log.info("Completed " + sample.getName());
	}

	private static final String[] FACEAWAY_HEADERS = {"Batch","Sample","Region","Size","Chr","ReadName",
            "RevStart","RevEnd","RevSeq","RevCigar","RevQuality",
            "FwdStart","FwdEnd","FwdSeq","FwdCigar","FwdQuality"};

	public class SampleResult  {
		Sample     sample;
		Faceaway[] faceaways;
		
		public SampleResult(Sample sample, Faceaway[] faceaways) {
			this.sample = sample;
			this.faceaways = faceaways;
		}
		
		/*
		 * Write out the results for this sample.
		 */
		protected void writetoFile (File outFolder) throws AnalysisException, IOException  {
			Sample sample = this.sample;

			TableOutput out = new TableOutput (outFolder, sample.getName()+".faceaway.tab", FACEAWAY_HEADERS, 64 * 1024);		
			Faceaway[] ff = this.faceaways;
			for (int fIdx = 0; fIdx < ff.length; fIdx++) {
				Faceaway f = ff[fIdx];
				out.newRow();
				out.appendValue(sample.getBatch());
				out.appendValue(sample.getName());
				out.appendValue(f.getRegion().toString());
				out.appendValue(f.getSize());
				out.appendValue(f.getRegion().getChromosome());
				out.appendValue(f.reverseRead.getReadName());
				
				out.appendValue(f.reverseRead.getAlignmentStart());
				out.appendValue(f.reverseRead.getAlignmentEnd());
				out.appendValue(f.reverseRead.getReadString());
				out.appendValue(f.reverseRead.getCigarString());
				out.appendValue(f.reverseRead.getMappingQuality());
				
				out.appendValue(f.forwardRead.getAlignmentStart());
				out.appendValue(f.forwardRead.getAlignmentEnd());
				out.appendValue(f.forwardRead.getReadString());
				out.appendValue(f.forwardRead.getCigarString());
				out.appendValue(f.forwardRead.getMappingQuality());
			}
			out.close();
		}
	}

	/* **********************************************************************
	 * Collective analysis of results from all samples
	 * **********************************************************************
	 */
	protected void analyzeAllSampleResults (Sample[] samples) throws AnalysisException, IOException  {
		mergeResultFiles(samples, "faceaway");
	}

	private void mergeResultFiles(Sample[] samples, String filenameSuffix) throws AnalysisException, IOException {
		TableOutput mergeOut = null;
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			File sampleFolder = getSampleSubfolder(outRootFolder, sample, true);
			File sampleFile = new File(sampleFolder, String.valueOf(sample.getName()) + "." + filenameSuffix + ".tab");
			if (!sampleFile.exists() || !sampleFile.canRead()) {
				log.warn("Could not access file " + sampleFile.getAbsolutePath() + " - skipping sample.");
			} else {
				TableInput tif = new TableInput(sampleFile);
				// Get the headers, removing the "Num" automatic field at the start
				String[] fieldHeaders = tif.getFieldNames();
				fieldHeaders = Arrays.copyOfRange(fieldHeaders, 1, fieldHeaders.length); 
				if (mergeOut == null) {
					mergeOut = new TableOutput(outRootFolder, "AllSamples." + filenameSuffix + ".tab", fieldHeaders, 64 * 1024);
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
			if (args.length < 4) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$SingleSample <batchId> <sampleId> <bamFile> <rootFolder>");
				log.error("       if <batchId> is \"NULL\", use 4-letter file prefix as batch");
				log.error("       if <batchId> is \"-\", use NO_BATCH as group for this sample");
				return;
			}
			String batchId = "NULL".equals(args[0]) ? null : args[0];	log.info("BatchId: "+batchId);
			String sampleId = args[1];									log.info("SampleId: "+sampleId);
			File sampleBamFile = new File(args[2]);	    				log.info("SampleBamFile: "+sampleBamFile.getAbsolutePath());
			File rootFolder = new File(args[3]);						log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				FaceawayAnalysis task = new FaceawayAnalysis(rootFolder);
				Sample sample = new Sample (batchId, sampleId, sampleBamFile);
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
			if (args.length < 2) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MultiSample <sampleListFile> <rootFolder>");
				return;
			}
			File sampleListFile = new File(args[0]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File rootFolder = new File(args[1]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());

			int maxThreads = Integer.parseInt(System.getProperty("maxThreads","0"));
			
			try {	
				FaceawayAnalysis task = new FaceawayAnalysis(rootFolder);
				MultiSampleAnalysis multi = new MultiSampleAnalysis(sampleListFile, maxThreads);
				multi.execute(task);
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
			if (args.length < 2) {
				log.error("Usage: org.cggh.bam.sampleClass.SampleClassAnalysis$MergeResults <sampleListFile> <rootFolder>");
				return;
			}
			File sampleListFile = new File(args[0]);	log.info("SampleListFile: "+sampleListFile.getAbsolutePath());
			File rootFolder = new File(args[1]);		log.info("RootFolder: "+rootFolder.getAbsolutePath());
			
			try {
				Sample[] samples = new SampleList(sampleListFile, false).getSamples();
				FaceawayAnalysis task = new FaceawayAnalysis(rootFolder);
				task.analyzeAllSampleResults(samples);
			} catch (Exception e) {
				log.error("Error executing task: " + e);
				return;
			}
			log.info("Exiting");
		}
	}
}
