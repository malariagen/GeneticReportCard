package org.cggh.bam.crtGenotyping.mk1;

import org.cggh.common.counters.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.sequence.*;
import org.cggh.common.util.*;
import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;


public class CrtHaploAnalysis {
	
	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private File samFolder;
	private File outFolder;
	
	private AnchorSeq[] crtAnchors = new AnchorSeq[] {
			new AnchorSeq (403593, "TATTATTTATTTAAGTGTA"),
			new AnchorSeq (403627, "ATTTTTGCTAAAAGAAC"),
	};
	private String WILD_TYPE_HAPLO = "CVMNK";
	
	private int haploStart = 403612;
	private int haploEnd   = 403626;
	
	private int bufferSize =  4 * 1024 * 1024;
	
	/* ==========================================================
	 * Samples and SNP lists initialization
	 * ==========================================================
	 */
	public CrtHaploAnalysis (File samFolder, File outFolder) throws AnalysisException  {
		this.samFolder = samFolder;
		this.outFolder = outFolder;
		if (!outFolder.exists()) {
			outFolder.mkdirs();
		}
	}

	HaploCounters haploSampleCounters = new HaploCounters();
	
	public void execute () throws AnalysisException, IOException  {
		log.info("SAM files folder: "+samFolder);
		final String fileSuffix = "-crt-flankReads.sam";		
		String[] fnames = samFolder.list(new FilenameFilter() {
			    public boolean accept(File dir, String name) {
			    	return name.endsWith(fileSuffix);
			    }});
		
		ArrayList<SampleSummary> sampleSummaryList = new ArrayList<SampleSummary>();
		for (int sIdx = 0; sIdx < fnames.length; sIdx++) {
			// Get sample name and file
			String sampleName = fnames[sIdx].substring(0, fnames[sIdx].length()-fileSuffix.length());
			File samFile = new File (samFolder, fnames[sIdx]);
			
			// Read the reads from the SAM file
			SeqRead[] sampleReads = getSampleReads (samFile);
			
			// Perform a reads alignment and count the haplotypes in the specified region
			ReadAlignment ra = new ReadAlignment(sampleName, sampleReads);
			log.info(ra.name+": " + ra.sampleReads.length + " reads; alignment: "+ra.alignStart+"-"+(ra.alignStart+ra.alignLen-1));		
			sampleSummaryList.add(ra.getSampleSummary());

			// Write out to file the reads alignment 
			File fastaFile = new File(outFolder, sampleName+".fasta");
			ra.writeFastaFile(fastaFile);
			
			// Record which haplotypes we found in a master list
			LabelCounter[] hc = ra.haploCounters.getSortedCounters();
			for (int hIdx = 0; hIdx < hc.length; hIdx++) {
				haploSampleCounters.increment(hc[hIdx].getLabel());
			}
		}
		SampleSummary[] sampleSummaries = sampleSummaryList.toArray(new SampleSummary[sampleSummaryList.size()]);
		outputSampleSummaries (sampleSummaries);

		// Compute some haplotype statistics
		HaplotypeStats[] hStats = computeHaplotypeStats (haploSampleCounters, sampleSummaries);
		
		// Write out a haplotype read count table
		outputHaploReadCounts (hStats, sampleSummaries, "-beforeMinAlt");
		
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
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			ss.haploCounters.filterCounters(haploFilter);
		}
		
		// Final step: make a call only if we have at least 5 reads for a valid haplo
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			int readCount = ss.haploCounters.getTotal();
			if (readCount < 5) {
				ss.haploCounters.clear();
			}
		}
		
		// Write out a final haplotype read count table
		outputHaploReadCounts (filteredHStats, sampleSummaries, "-final");
		
		// Write out a haplotype stats table
		outputHaploStats (filteredHStats, "-final");
		
		// Write out the final calls
		 outputSampleHaplotypeCalls (sampleSummaries);
	}
	
	private HaplotypeStats[] computeHaplotypeStats (HaploCounters haploSampleCounters, SampleSummary[] sampleSummaries) {
		LabelCounter[] hSampleCounters = haploSampleCounters.getSortedCounters();
		HaplotypeStats[] hStats = new HaplotypeStats[hSampleCounters.length];
		for (int hIdx = 0; hIdx < hSampleCounters.length; hIdx++) {
			HaploCounter hc = (HaploCounter)hSampleCounters[hIdx];
			hStats[hIdx] = new HaplotypeStats(hc.getLabel(), hc.fullLabel);
		}
		
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			double totalSampleReads = ss.haploCounters.getTotal();
			for (int hIdx = 0; hIdx < hSampleCounters.length; hIdx++) {
				HaplotypeStats stat = hStats[hIdx];
				LabelCounter c = ss.haploCounters.getCounter(stat.haplo);
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
	
	
	private SeqRead[] getSampleReads (File samFile) throws AnalysisException {

		ArrayList<SeqRead> readsList = new ArrayList<SeqRead>();
		TabFileReader tabReader = null;
		try {
			tabReader = new TabFileReader(new FileReader(samFile));
			String[] fields;
			while (true) {
				fields = tabReader.getNextValidLine();
				if (fields == null) {
					break;
				}
				SeqRead sr = new SeqRead(fields);
				readsList.add(sr);
			}
		} catch (IOException e) {
			throw new AnalysisException ("Error reading file "+ samFile.getAbsolutePath()+ ": "+ e);
		} finally {
			try {
				tabReader.close();
			} catch (IOException e) {
				throw new AnalysisException ("Error closing file "+ samFile.getAbsolutePath()+ ": "+ e);
			}
		}
		return readsList.toArray(new SeqRead[readsList.size()]);
	}

	public class ReadAlignment {
		String    name;
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
		
		public ReadAlignment (String name, SeqRead[] sampleReads) {
			this.name = name;
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
		
		public SampleSummary getSampleSummary () {
			return new SampleSummary (name, sampleReads.length, haploCount, misalignCount, rawHaploCounters, haploCounters);
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
	
	public class SampleSummary {
		String    name;
		int       readCount;
		int       haploCount;
		int       misalignCount;
		LabelCounters rawHaploCounters;
		LabelCounters haploCounters;
		
		public SampleSummary(String name, int readCount, int haploCount, int misalignCount, LabelCounters rawHaploCounters, LabelCounters haploCounters) {
			this.name = name;
			this.readCount = readCount;
			this.haploCount = haploCount;
			this.misalignCount = misalignCount;
			this.rawHaploCounters = rawHaploCounters;
			this.haploCounters = haploCounters;
		}
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
			this.id = samFields[0];
			this.sequence = samFields[9];
			this.quality = samFields[10];
			
			anchorSequence(sequence);
			if (anchoredSequence == null) {
				anchorSequence(getReverseComplementSeq());
			}
			
			if (anchoredSequence != null) {
				int haploLen = 1 + haploEnd - haploStart;
				int readLen = anchoredSequence.length();
				int haploStartOffset = haploStart - startPos;
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
		
		public String  getReverseComplementSeq() {
			int len = sequence.length();
			char[] result = new char[len];
			for (int i = 0; i < len; i++) {
				char nt = sequence.charAt(len-i-1);
				result[i] = ((nt == 'A') ? 'T' : ((nt == 'T') ? 'A' : ((nt == 'C') ? 'G' : 'C')));
			}
			return new String(result);
		}
	}
	 
	public class AnchorSeq {
		int pos;
		String sequence;
		
		public AnchorSeq(int pos, String sequence) {
			this.pos = pos;
			this.sequence = sequence;
		}
	}
	 
	/* ==========================================================
	 * 
	 * ==========================================================
	 */
	private void outputSampleSummaries (SampleSummary[] sampleSummaries) throws AnalysisException {
		// Open the sample summary file
		String[] colNames = {"Sample","Haplotypes","UnfilteredHaplotypes","TotalReads","HaplotypeCount","MisalignedReads"};
		TableOutput sampleSummaryOut = new TableOutput (outFolder, "HaplotypeSummariesBySample-beforeMinAlt.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			sampleSummaryOut.newRow();
			sampleSummaryOut.appendValue(ss.name);
			sampleSummaryOut.appendValue(ss.haploCounters.getSummary());
			sampleSummaryOut.appendValue(ss.rawHaploCounters.getSummary());
			sampleSummaryOut.appendValue(ss.readCount);
			sampleSummaryOut.appendValue(ss.haploCount);
			sampleSummaryOut.appendValue(ss.misalignCount);
		}
		sampleSummaryOut.close();
	}

	private void outputSampleHaplotypeCalls (SampleSummary[] sampleSummaries) throws AnalysisException {
		// Open the sample summary file
		String[] colNames = {"Sample","Call","Haplotypes","HaplotypeReads"};
		
		TableOutput sampleSummaryOut = new TableOutput (outFolder, "HaplotypeCallsBySample-final.tab", colNames, bufferSize);
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			sampleSummaryOut.newRow();
			sampleSummaryOut.appendValue(ss.name);
			
			if (ss.haploCounters.isEmpty()) {
				sampleSummaryOut.appendValue("MI");
				sampleSummaryOut.appendBlankValue();
				sampleSummaryOut.appendBlankValue();
			} else {
				LabelCounter[] counters = ss.haploCounters.getSortedCounters();
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
			sampleSummaryOut.appendValue(ss.haploCounters.getSummary());
		}
		sampleSummaryOut.close();
	}

	private void outputHaploReadCounts (HaplotypeStats[] hStats, SampleSummary[] sampleSummaries, String suffix) throws AnalysisException {
		// Write out a haplotype read count table
		String[] hDisplayLabels = new String[hStats.length];
		for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
			hDisplayLabels[hIdx] = hStats[hIdx].displayLabel;
		}
		TableOutput haploReadsOut = new TableOutput (outFolder, "HaplotypeSampleCount"+suffix+".tab", 
				               TextUtilities.mergeStringLists(new String[] {"Sample"}, hDisplayLabels), bufferSize);
		for (int sIdx = 0; sIdx < sampleSummaries.length; sIdx++) {
			SampleSummary ss = sampleSummaries[sIdx];
			haploReadsOut.newRow();
			haploReadsOut.appendValue(ss.name);
			for (int hIdx = 0; hIdx < hStats.length; hIdx++) {
				LabelCounter c = ss.haploCounters.getCounter(hStats[hIdx].haplo);
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
			System.out.println("Usage: org.cggh.bam.genotyping.crt.mk1.CrtHaploAnalysisPrev <samFolder> <outFolder>");
			return;
		}
		File samFolder = new File(args[0]);
		File outFolder = new File(args[1]);
			
		try {
			CrtHaploAnalysis task = new CrtHaploAnalysis(samFolder, outFolder);
			task.execute();
		} catch (Exception e) {
			log.error("Error: " + e);
			return;
		}
		log.info("Exiting");
	}

}
