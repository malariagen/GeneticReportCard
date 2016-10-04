package org.cggh.bam.genotyping.adhoc;

import org.cggh.bam.genotyping.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.cggh.common.sequence.io.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class AlignmentRefFixer {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private Sample[]      samples;
	private Locus[] locusConfigs;
	
	private File lociFolder;
	
	
	/* ==========================================================
	 * Samples and SNP lists initialization
	 * ==========================================================
	 */
	public AlignmentRefFixer (File configFile, File sampleListFile, File refFastaFile, File rootFolder) throws AnalysisException  {

		// Load up the reference genome sequences
		ReferenceGenome.initialize(refFastaFile);

		GenotypingConfig config = new GenotypingConfig (configFile);
		locusConfigs = config.getLoci();
		
		// Get the list of samples
		samples = readSamples(sampleListFile);

		// Create the root output directory
		this.lociFolder    = FileUtilities.checkFolder(rootFolder, "loci", false);
	}


	public void execute () throws AnalysisException, IOException  {
		
		for (int lIdx = 0; lIdx < locusConfigs.length; lIdx++) {
			
			log.info("Analyzing locus "+locusConfigs[lIdx].getName());
			ReadsAlignmentRetriever raRetriever = new ReadsAlignmentRetriever (locusConfigs[lIdx], lociFolder);

			for (int sIdx = 0; sIdx < samples.length; sIdx++) {
				ReadsAlignment ra = raRetriever.retrieveReadsAlignment(samples[sIdx]);
				if (ra == null) {
					continue;
				}
				int raStartPos = ra.getAlignStart();
				if (raStartPos < 0) {
					continue;
				}
				
				File fastaFile = raRetriever.getFastaFile(samples[sIdx]);
				Sequence[] raSeqs = null;
				try {
					SequenceSetReader seqSetReader = new SequenceSetReader (new FastaSequenceReader());
					raSeqs = seqSetReader.readSequences(fastaFile);
				} catch (Exception e) {
					log.error ("Error reading reads alignment FASTA file for sample "+samples[sIdx].getName()+": "+e);
					continue;
				}
				if (raSeqs.length == 0) {
					// Empty alignment
					continue;
				}
				if (raSeqs[0].getId().startsWith("REF|")) {
					// Already has a reference
					continue;
				}

				int raStartLen = ra.getAlignLen();
				String chrName = locusConfigs[lIdx].getReadSearchInterval().getChromosome();
				Sequence chrSeq = ReferenceGenome.getChrSequence(chrName);
				String refSeqData = chrSeq.getData().substring(raStartPos-1, raStartPos+raStartLen-1);
				String refSeqTitle = "REF|"+chrName+":"+raStartPos+"-"+(raStartPos+raStartLen-1);
				
				Sequence[] newSeqs = new Sequence[raSeqs.length+1];
				newSeqs[0] = new Sequence(refSeqTitle, refSeqData);
				System.arraycopy(raSeqs, 0, newSeqs, 1, raSeqs.length);
				
				// Write out to file the reads alignment
				String fasta = SequenceUtilities.makeFastaAlignment(newSeqs);
				try {
					FileUtilities.writeFileContent(fasta, fastaFile);
				} catch (IOException e) {
					throw new AnalysisException ("Error writing file "+ fastaFile.getAbsolutePath()+ ": "+ e);
				}
			}
		}
	}	
	
	/* *************************************************************************
	 * Initialization routines
	 * *************************************************************************
	 */
	private Sample[] readSamples(File sampleListFile) throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(sampleListFile));
		String[] colNames = new String[]{"Sample"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample(values[0], null);
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
			System.out.println("Usage: org.cggh.bam.genotyping.adhoc.AlignmentRefFixer <configFile> <sampleListFile> <refFasta> <rootFolder>");
			return;
		}
		File configFile = new File(args[0]);
		File sampleListFile = new File(args[1]);
		File refFastaFile = new File(args[2]);
		File rootFolder = new File(args[3]);
		
		try {
			AlignmentRefFixer task = new AlignmentRefFixer(configFile, sampleListFile, refFastaFile, rootFolder);
			task.execute();
		} catch (Exception e) {
			log.error("Error: " + e);
			return;
		}
		log.info("Exiting");
	}

}
