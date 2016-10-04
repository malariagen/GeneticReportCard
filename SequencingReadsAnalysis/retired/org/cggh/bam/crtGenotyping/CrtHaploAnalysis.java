package org.cggh.bam.crtGenotyping;

import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.GenomeRegion;
import org.cggh.common.sequence.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;

import htsjdk.samtools.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;


public class CrtHaploAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private File bamListFile;
	private File outFolder;
	
	private AnchorSeq[] crtAnchors = new AnchorSeq[] {
			new AnchorSeq (403593, "TATTATTTATTTAAGTGTA"),
			new AnchorSeq (403627, "ATTTTTGCTAAAAGAAC"),
	};
	private String WILD_TYPE_HAPLO = "CVMNK";
	
	private GenomeRegion haploSpan = GenomeRegion.parseRegion("Pf3D7_07_v3:403612-403626");
	
	private GenomeRegion readSearchInterval = GenomeRegion.parseRegion("Pf3D7_07_v3:403000-404000");
	
	private int bufferSize =  4 * 1024 * 1024;
	
	/* ==========================================================
	 * Samples and SNP lists initialization
	 * ==========================================================
	 */
	public CrtHaploAnalysis (File bamListFile, File outFolder) throws AnalysisException  {
		this.bamListFile = bamListFile;
		this.outFolder = outFolder;
		if (!outFolder.exists()) {
			outFolder.mkdirs();
		}
	}

	
	HaploCounters haploSampleCounters = new HaploCounters();
	
	public void execute () throws AnalysisException, IOException  {
		
		Sample[] samples = readSamples(bamListFile);
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {

			// Read the reads from the SAM file
			Sample sample = samples[sIdx];
			SeqRead[] sampleReads = getSampleReads (sample);
			
			// Perform a reads alignment and count the haplotypes in the specified region
			ReadAlignment ra = new ReadAlignment(sampleReads);
			sample.readCount = ra.sampleReads.length;
			sample.haploCount = ra.haploCount;
			sample.misalignCount = ra.misalignCount;
			sample.rawHaploCounters = ra.rawHaploCounters;
			sample.haploCounters = ra.haploCounters;
			log.info(sample.name+": " + sample.readCount + " reads; alignment: "+ra.alignStart+"-"+(ra.alignStart+ra.alignLen-1));

			// Write out to file the reads alignment 
			File fastaFile = new File(outFolder, sample.name+".fasta");
			ra.writeFastaFile(fastaFile);
			
			// Record which haplotypes we found in a master list
			LabelCounter[] hc = ra.haploCounters.getSortedCounters();
			for (int hIdx = 0; hIdx < hc.length; hIdx++) {
				haploSampleCounters.increment(hc[hIdx].getLabel());
			}
		}
		outputSampleSummaries (samples);

		// Compute some haplotype statistics
		HaplotypeStats[] hStats = computeHaplotypeStats (haploSampleCounters, samples);
		
		// Write out a haplotype read count table
		outputHaploReadCounts (hStats, samples, "-beforeMinAlt");
		
		// Write out a haplotype stats table
		outputHaploStats (hStats, "-beforeMinAlt");
		
		// Now go through and filter the haplotypes to reject likely artefact
		// The test is MinAlt: homozygous or >=10 reads in at least one sample
		ArrayList<HaplotypeStats> filteredHStatList = new ArrayList<HaplotypeStats>();
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			if ((hStats[hIdx].maxReadsFraction == 1.0) || (hStats[hIdx].maxReads >= 10)) {
				filteredHStatList.add(hStats[hIdx]);
			}
		}
		HaplotypeStats[] filteredHStats = filteredHStatList.toArray(new HaplotypeStats[filteredHStatList.size()]);
		String[] filteredHaplos = new String[filteredHStats.length];
		for (int hIdx = 0; hIdx < filteredHStats.length; hIdx++) {
			filteredHaplos[hIdx] = filteredHStats[hIdx].haplo;
		}
		
		// Go through the samples, keeping only the haplotypes that passed the MinAltfilter.
		LabelCounterFilter haploFilter = new LabelCounterFilterByLabelList (filteredHaplos);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			samples[sIdx].haploCounters.filterCounters(haploFilter);
		}
		
		// Final step: make a call only if we have at least 5 reads for a valid haplo
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			int readCount = samples[sIdx].haploCounters.getTotal();
			if (readCount < 5) {
				samples[sIdx].haploCounters.clear();
			}
		}
		
		// Write out a final haplotype read count table
		outputHaploReadCounts (filteredHStats, samples, "-final");
		
		// Write out a haplotype stats table
		outputHaploStats (filteredHStats, "-final");
		
		// Write out the final calls
		 outputSampleHaplotypeCalls (samples);
	}
	
	private Sample[] readSamples(File bamListFile) throws AnalysisException {
		ArrayList<Sample> sampleList = new ArrayList<Sample>();
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(bamListFile));
		String[] colNames = new String[]{"Sample","BamFile"};
		ColumnFileReader.ColumnReader cr = cfr.getColumnReader(colNames);
		while (cfr.nextRecord()) {
			String[] values = cr.getValues();
			Sample s = new Sample();
			s.name = values[0];
			s.bamFile = new File(values[1]);
			if (!s.bamFile.canRead()) {
				throw new AnalysisException("Error getting samples list: BAM file " + values[1] + " cannot be read.");
			}
			sampleList.add(s);
		}
		cfr.close();
		return sampleList.toArray(new Sample[sampleList.size()]);
	}

	
	private HaplotypeStats[] computeHaplotypeStats (HaploCounters haploSampleCounters, Sample[] samples) {
		LabelCounter[] hSampleCounters = haploSampleCounters.getSortedCounters();
		HaplotypeStats[] hStats = new HaplotypeStats[hSampleCounters.length];
		for (int hIdx = 0; hIdx < hSampleCounters.length; hIdx++) {
			HaploCounter hc = (HaploCounter)hSampleCounters[hIdx];
			hStats[hIdx] = new HaplotypeStats(hc.getLabel(), hc.fullLabel);
		}
		
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			double totalSampleReads = samples[sIdx].haploCounters.getTotal();
			for (int hIdx = 0; hIdx < hSampleCounters.length; hIdx++) {
				HaplotypeStats stat = hStats[hIdx];
				LabelCounter c = samples[sIdx].haploCounters.getCounter(stat.haplo);
				if (c != null) {
					int readCount = c.getCount();
					double readsFraction = ((double)readCount) / totalSampleReads;
					stat.sampleCount++;
					stat.maxReads = (readCount > stat.maxReads) ? readCount : stat.maxReads;
					stat.maxReadsFraction = (readsFraction > stat.maxReadsFraction) ? readsFraction : stat.maxReadsFraction;
				}
			}
		}
		return hStats;
	}
	
	
	public static class LabelCounterFilterByLabelList implements LabelCounterFilter {
		private String[] labelsAllowed;
		public LabelCounterFilterByLabelList (String[] labelsAllowed) {
			this.labelsAllowed = labelsAllowed;
		}
		public boolean isCounterValid (LabelCounter counter) {
			String l = counter.getLabel();
			for (int i = 0; i < labelsAllowed.length; i++) {
				if (labelsAllowed[i].equals(l)) {
					return true;
				}
			}
			return false;
		}
	}
	
	public class HaplotypeStats {
		String haplo;
		String displayLabel;
		int    sampleCount;
		int    maxReads;
		double maxReadsFraction;
		
		public HaplotypeStats(String haplo, String displayLabel) {
			super();
			this.haplo = haplo;
			this.displayLabel = displayLabel;
		}
	}
	
	
	SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	private SeqRead[] getSampleReads (Sample sample) throws AnalysisException {
		
		// Build the regexes that will be used to identify the useful reads
		String anchorRegex = crtAnchors[0].sequence    + "|" + crtAnchors[1].sequence 
				     + "|" + crtAnchors[0].revSequence + "|" + crtAnchors[1].revSequence;
	    Pattern anchorPattern = Pattern.compile(anchorRegex);
		
		ArrayList<SeqRead> readsList = new ArrayList<SeqRead>();
		SamReader samReader = samReaderFactory.open(sample.bamFile);
		
		// First search reads aligned onto the pfcrt region
		filterSampleReads (samReader.query(readSearchInterval.getChromosome(), readSearchInterval.getStartPos(), readSearchInterval.getStopPos(), true), 
				anchorPattern, readsList);

		// Then search unmapped reads (this does not necessarily find them all, but should do)
		filterSampleReads (samReader.queryUnmapped(), anchorPattern, readsList);

		return readsList.toArray(new SeqRead[readsList.size()]);
	}
	
	private void filterSampleReads (SAMRecordIterator it, Pattern regexPattern, ArrayList<SeqRead> readsList) throws AnalysisException {
		while (it.hasNext()) {
			SAMRecord record = it.next();
			String readSeq = record.getReadString();
		    if (regexPattern.matcher(readSeq).find()) {
				SeqRead sr = new SeqRead(record.getReadName(), readSeq, record.getBaseQualityString());
				readsList.add(sr);
		    }
		}
		it.close();
	}
	
	
	public class ReadAlignment {
		SeqRead[] sampleReads;
		int       alignStart;
		int       alignLen;
		char[][]  alignment;
		char[]    consensus;
		int[]     differences;
		int       haploCount;
		int       misalignCount;
		HaploCounters rawHaploCounters;
		HaploCounters haploCounters;
		
		public ReadAlignment (SeqRead[] sampleReads) {
			this.sampleReads = sampleReads;
			
			// Align all reads according to their anchor position
			alignment = makeAlignment ();
			
			// Get the consensus sequence for this alignment
			consensus = computeConsensus();
			
			// Determine the count of haplotypes in this alignment
			haploCount = 0;
			misalignCount = 0;
			rawHaploCounters = new HaploCounters();
			haploCounters = new HaploCounters();
			
			for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
				SeqRead sr = sampleReads[rIdx];
			
				// Determine the number of differences wrt the consensus
				sr.differences = 0;
				for (int i = 0; i < alignLen; i++) {
					char nt = alignment[rIdx][i];
					if ((nt != '-') && (nt != consensus[i])) {
						sr.differences++;
					}
				}

				// If there are too many, we discard the sequence as misaligned
				if (sr.hasTooManyDifferences()) {
					misalignCount++;
					continue;
				}
				
				// Does this read contain the full haplotype? If so, catalogue it.
				String haplo = sr.ntHaplotype;
				if (haplo != null) {
					haploCount++;
					rawHaploCounters.increment(haplo);
					if (!sr.hasLowQualityHaplo()) {
						haploCounters.increment(haplo);
					}
				}
			}
			
			LabelCounterFilter nRemover = new LabelCounterFilterByLabelSubstring("N");
			haploCounters.filterCounters(nRemover);
			
			LabelCounterFilter singletonRemover = new LabelCounterFilterByMinCount(2);
			haploCounters.filterCounters(singletonRemover);
		}
				
		private char[][] makeAlignment () {
			// Compute the maximum extent of the alignemnt
			int alignEnd = alignStart = -1;
			for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
				SeqRead sr = sampleReads[rIdx];
				int rStart = sr.startPos;
				int rEnd = sr.startPos + sr.sequence.length() - 1;
				if (rIdx == 0) {
					alignStart = rStart;
					alignEnd = rEnd;
				} else {
					alignStart = (alignStart < rStart) ? alignStart : rStart;
					alignEnd = (alignEnd > rEnd) ? alignEnd : rEnd;
				}
			}
			alignLen = 1 + alignEnd - alignStart;

			// Create an array for storing the alignments
			char[][] result = new char[sampleReads.length][alignLen];
			
			// Write the read sequences in the right places
			for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
				Arrays.fill(result[rIdx], '-');
				SeqRead sr = sampleReads[rIdx];
				int idx = sr.startPos - alignStart;
				for (int i = 0; i < sr.anchoredSequence.length(); i++) {
					result[rIdx][idx+i]=sr.anchoredSequence.charAt(i);
				}
			}
			return result;
		}
		
		private char[] computeConsensus () {
			char[] result = new char[alignLen];
			int[] ntCounts = new int[4];
			for (int i = 0; i < alignLen; i++) {
				ntCounts[0] = ntCounts[1] = ntCounts[2] = ntCounts[3] = 0;
				for (int j = 0; j < sampleReads.length; j++) {
					char nt = alignment[j][i];
					if ((nt != '-') && (nt != 'N')) {
						int ntIdx = nt2Index (nt);
						ntCounts[ntIdx]++;
					}
				}
				result[i] = getMajorityAllele (ntCounts);
			}
			return result;
		}
		
		private char getMajorityAllele (int[] ntCounts) {
			return (index2Nt(getMajorityIndex(ntCounts)));
		}
		
		private int getMajorityIndex (int[] ntCounts) {
			int mCount = -1;
			int mIdx = -1;
			for (int i = 0; i < 4; i++) {
				if (ntCounts[i] > mCount) {
					mIdx = i;
					mCount = ntCounts[i];
				}
			}
			return mIdx;
		}
		
		private int nt2Index (char nt) {
			switch (nt) {
			case 'A': return 0;
			case 'C': return 1;
			case 'G': return 2;
			case 'T': return 3;
			}
			return -1;
		}
		private char index2Nt (int index) {
			switch (index) {
			case 0: return 'A';
			case 1: return 'C';
			case 2: return 'G';
			case 3: return 'T';
			}
			return (char)0;
		}
		
		public void writeFastaFile (File fastaFile) throws AnalysisException {
			Sequence[] alignSeq = new Sequence[sampleReads.length];
			for (int rIdx = 0; rIdx < sampleReads.length; rIdx++) {
				SeqRead sr = sampleReads[rIdx];
				String seqTitle = sr.id;
				if (sr.hasTooManyDifferences()) {
					seqTitle = "[MISALIGNED] " + seqTitle;
				}
				if (sr.ntHaplotype == null) {
					seqTitle = "[NO_HAPLO] " + seqTitle;
				}
				if (sr.hasLowQualityHaplo()) {
					seqTitle = "[LOW_HAPLO_Q] " + seqTitle;
				}
				alignSeq[rIdx] = new Sequence (seqTitle, new String(alignment[rIdx]));
			}
			String fasta = SequenceUtilities.makeFastaAlignment(alignSeq);
			try {
				FileUtilities.writeFileContent(fasta, fastaFile);
			} catch (IOException e) {
				throw new AnalysisException ("Error writing file "+ fastaFile.getAbsolutePath()+ ": "+ e);
			}
		}
	}
	
	private class Sample {
		String name;
		File bamFile;
		int       readCount;
		int       haploCount;
		int       misalignCount;
		LabelCounters rawHaploCounters;
		LabelCounters haploCounters;
	}
	
	/* ==========================================================
	 * 
	 * ==========================================================
	 */
	static final int MAX_DIFFERENCES_FROM_CONSENSUS = 10;
	static final int MIN_PHRED_SCORE = 20;
	
	public class SeqRead {
		String id;
		String sequence;
		String quality;
		String anchoredSequence;
		String ntHaplotype;
		int minHaploQuality;
		AnchorSeq anchor;
		int startPos;
		int anchorPos;
		int differences;
		
		public SeqRead(String[] samFields) {
			this (samFields[0], samFields[9], samFields[10]);
		}
		
		public SeqRead(String readId, String readString, String qualityString) {
			this.id = readId;
			this.sequence = readString;
			this.quality = qualityString;
			
			anchorSequence(sequence);
			if (anchoredSequence == null) {
				anchorSequence(getReverseComplementSeq(sequence));
			}
			
			if (anchoredSequence != null) {
				int haploLen = 1 + haploSpan.getStopPos() - haploSpan.getStartPos();
				int readLen = anchoredSequence.length();
				int haploStartOffset = haploSpan.getStartPos() - startPos;
				int haploEndOffset = haploStartOffset + haploLen - 1;
				if ((haploStartOffset >= 0) && (haploEndOffset < readLen)) {
					// Read covers the complete haplotype, extract it
					String h = anchoredSequence.substring(haploStartOffset, haploStartOffset+haploLen);
					if (!h.contains("N")) {
						ntHaplotype = h;
						minHaploQuality = 1000;
						for (int i = 0; i < haploLen; i++) {
							int q = getPhredScore(haploStartOffset+i);
							if (q < minHaploQuality) {
								minHaploQuality = q;
							}
						}
					}
				}
			}
		}
		
		public int getPhredScore (int offset) {
			int phredCode = quality.charAt(offset);
			return (phredCode - 33);
		}
		
		public boolean hasTooManyDifferences () {
			return (differences > MAX_DIFFERENCES_FROM_CONSENSUS);
		}
		
		public boolean hasLowQualityHaplo () {
			return (minHaploQuality < MIN_PHRED_SCORE);
		}
		
		private void anchorSequence (String seq) {
			for (int i = 0; i < crtAnchors.length; i++) {
				AnchorSeq a = crtAnchors[i];
				int pos = seq.indexOf(a.sequence);
				if (pos >= 0) {
					anchoredSequence = seq;
					anchor = a;
					anchorPos = pos;
					startPos = a.pos - pos;
					break;
				}
			}
		}
	}
	 
	public class AnchorSeq {
		int pos;
		String sequence;
		String revSequence;
		
		public AnchorSeq(int pos, String sequence) {
			this.pos = pos;
			this.sequence = sequence;
			this.revSequence = getReverseComplementSeq(sequence);
		}
	}
	
	public static String  getReverseComplementSeq(String sequence) {
		int len = sequence.length();
		char[] result = new char[len];
		for (int i = 0; i < len; i++) {
			char nt = sequence.charAt(len-i-1);
			result[i] = ((nt == 'A') ? 'T' : ((nt == 'T') ? 'A' : ((nt == 'C') ? 'G' : 'C')));
		}
		return new String(result);
	}
	 
	/* ==========================================================
	 * 
	 * ==========================================================
	 */
	private void outputSampleSummaries (Sample[] samples) throws AnalysisException {
		// Open the sample summary file
		String[] colNames = {"Sample","Haplotypes","UnfilteredHaplotypes","TotalReads","HaplotypeCount","MisalignedReads"};
		TableOutput sampleSummaryOut = new TableOutput (outFolder, "HaplotypeSummariesBySample-beforeMinAlt.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			sampleSummaryOut.newRow();
			sampleSummaryOut.appendValue(sample.name);
			sampleSummaryOut.appendValue(sample.haploCounters.getSummary());
			sampleSummaryOut.appendValue(sample.rawHaploCounters.getSummary());
			sampleSummaryOut.appendValue(sample.readCount);
			sampleSummaryOut.appendValue(sample.haploCount);
			sampleSummaryOut.appendValue(sample.misalignCount);
		}
		sampleSummaryOut.close();
	}

	private void outputSampleHaplotypeCalls (Sample[] samples) throws AnalysisException {
		// Open the sample summary file
		String[] colNames = {"Sample","Call","Haplotypes","HaplotypeReads"};
		
		TableOutput sampleSummaryOut = new TableOutput (outFolder, "HaplotypeCallsBySample-final.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			sampleSummaryOut.newRow();
			sampleSummaryOut.appendValue(sample.name);
			
			if (sample.haploCounters.isEmpty()) {
				sampleSummaryOut.appendValue("MI");
				sampleSummaryOut.appendBlankValue();
				sampleSummaryOut.appendBlankValue();
			} else {
				LabelCounter[] counters = sample.haploCounters.getSortedCounters();
				if (counters.length == 1) {
					String aaHaplo = ((HaploCounter)counters[0]).aaSequence;
					sampleSummaryOut.appendValue(WILD_TYPE_HAPLO.equals(aaHaplo) ? "WT" : "MU");
					sampleSummaryOut.appendValue(aaHaplo);
				} else {
					sampleSummaryOut.appendValue("<Het>");
					StringBuffer sb = new StringBuffer();
					for (int cIdx = 0; cIdx < counters.length; cIdx++) {
						if (cIdx > 0) {
							sb.append(",");
						}
						sb.append(((HaploCounter)counters[cIdx]).aaSequence);
					}
					sampleSummaryOut.appendValue(sb.toString());
				}
			}
			sampleSummaryOut.appendValue(sample.haploCounters.getSummary());
		}
		sampleSummaryOut.close();
	}

	private void outputHaploReadCounts (HaplotypeStats[] hStats, Sample[] samples, String suffix) throws AnalysisException {
		// Write out a haplotype read count table
		String[] hDisplayLabels = new String[hStats.length];
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			hDisplayLabels[hIdx] = hStats[hIdx].displayLabel;
		}
		TableOutput haploReadsOut = new TableOutput (outFolder, "HaplotypeSampleCount"+suffix+".tab", 
				               TextUtilities.mergeStringLists(new String[] {"Sample"}, hDisplayLabels), bufferSize);
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			Sample sample = samples[sIdx];
			haploReadsOut.newRow();
			haploReadsOut.appendValue(sample.name);
			for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
				LabelCounter c = sample.haploCounters.getCounter(hStats[hIdx].haplo);
				haploReadsOut.appendValue((c != null) ? c.getCount() : 0);
			}
		}
		haploReadsOut.close();
	}

	private void outputHaploStats (HaplotypeStats[] hStats, String suffix) throws AnalysisException {
		TableOutput haploStatsOut = new TableOutput (outFolder, "HaplotypeStats"+suffix+".tab", 
                new String[] {"Haplotype","SampleCount","MaxReads","MaxReadFraction"}, bufferSize);
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			haploStatsOut.newRow();
			haploStatsOut.appendValue(hStats[hIdx].displayLabel);
			haploStatsOut.appendValue(hStats[hIdx].sampleCount);
			haploStatsOut.appendValue(hStats[hIdx].maxReads);
			haploStatsOut.appendValue(hStats[hIdx].maxReadsFraction);
		}
		haploStatsOut.close();
	}


	/* ==========================================================
	 * 
	 * ==========================================================
	 */
	public static class HaploCounter extends LabelCounter {
		public String ntSequence;
		public String aaSequence;
		public String fullLabel;
		
		public HaploCounter (String ntSequence) {
			super (ntSequence);
			this.ntSequence = ntSequence;
			this.aaSequence = translateHaplotype(ntSequence);
			this.fullLabel = ntSequence+"["+aaSequence+"]";
		}
		
		public String getSummary() {
			return appendCountToLabel(fullLabel);
		}
		
		private String translateHaplotype (String haplo) {
			StringBuffer sb = new StringBuffer();
			int aaLen = haplo.length() / 3;
			for (int i = 0; i < aaLen; i++) {
				int offset = i * 3;
				Codon c = Translation.getCodon(haplo.charAt(offset), haplo.charAt(offset+1), haplo.charAt(offset+2));
				sb.append(c.getAmino());
			}
			return sb.toString();
		}
	}
	
	public static class HaploCounters extends LabelCounters {
		public LabelCounter createCounter (String label) {
			return new HaploCounter(label);
		}
	}
	

	/* ==========================================================
	 * Execution
	 * ==========================================================
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println("Usage: org.cggh.bam.genotyping.crt.CrtHaploAnalysis <bamListFile> <outFolder>");
			return;
		}
		File bamListFile = new File(args[0]);
		File outFolder = new File(args[1]);
			
		try {
			CrtHaploAnalysis task = new CrtHaploAnalysis(bamListFile, outFolder);
			task.execute();
		} catch (Exception e) {
			log.error("Error: " + e);
			return;
		}
		log.info("Exiting");
	}

}
