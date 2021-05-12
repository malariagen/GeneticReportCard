package org.cggh.bam.depth;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.genome.*;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;

import java.io.*;
import java.util.*;
import org.apache.commons.logging.*;

public class DepthAnalysis extends BaseAnalysis {
	
	private static Log log = LogFactory.getLog((String)ClassUtilities.getCurrentClassName());
	
	private static final int LOCUS_END_IGNORE_LEN = 120;
	
	private DepthLocus[]    loci;
	private DepthPosition[] dpos;
	
	public DepthAnalysis(File bedFile, File rootFolder) throws AnalysisException {
		super(rootFolder);
		loci = getLoci (bedFile);
		dpos = getDepthPositions (loci);
	}
	
	private TableOutput skipOut = null;
	private ArrayList<SampleDepth> sampleDepthList = new ArrayList<SampleDepth>();

	public void analyzeAllSamples(Sample[] samples) throws AnalysisException {
		// Read in all the data
		for (int sIdx = 0; sIdx < samples.length; sIdx++) {
			// Read read count data for all the positions
			SampleDepth sd = processSampleReadCounts (samples[sIdx]);
			// If some problem occurred, just record the sample and skip.
			if (sd == null) {
				continue;
			}
			sampleDepthList.add(sd);
			if (((sIdx+1) % 100) == 0) {
				log.info("Processed samples: "+(sIdx+1));
			}
		}
		if (skipOut != null) {
			skipOut.close();
		}
		SampleDepth[] depths = sampleDepthList.toArray(new SampleDepth[sampleDepthList.size()]);
		
		// Now write out the read count values
		log.info("Writing read count file");
		writeDepths(depths, "AllSamples.readCounts.tab", true);
		writeSampleStats(depths);
		
		// Normalize coverage
		log.info("Normalizing coverage");
		for (int dIdx = 0; dIdx < depths.length; dIdx++) {
			// Read read count data for all the positions
			SampleDepth sd = depths[dIdx];
			for (int pIdx = 0; pIdx < dpos.length; pIdx++) {
				if (sd.sampleMedian == 0.0) {
					sd.depth[pIdx] = Double.NaN;
				} else {
					sd.depth[pIdx] /= sd.sampleMedian;
				}
			}
		}
		
		// Now write out the normalized depth
		log.info("Writing normalized coverage file");
		writeDepths(depths, "AllSamples.normDepth.tab", false);
		
	}

	
	private SampleDepth processSampleReadCounts (Sample sample) throws AnalysisException {
			
		// Read depth data for all the positions
		SampleDepth sd = readSampleDepths(sample, dpos);
		
		// If some problem occurred, just record the sample and skip.
		if (sd == null) {
			if (skipOut == null) {
				skipOut = new TableOutput(outRootFolder, "SkippedSamples.tab", new String[] {"Sample"}, 4096);
			}
			skipOut.newRow();
			skipOut.appendValue(sample.getName());
			return null;
		}
		
		// Work out the medians
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			DepthLocus locus = loci[lIdx];
			int startPosIdx = locus.startPosIdx + LOCUS_END_IGNORE_LEN;
			int posCount = locus.length - (2 * LOCUS_END_IGNORE_LEN);
			
			Statistics dStat = new Statistics(sd.depth, startPosIdx, startPosIdx+posCount);
			sd.locusMedians[lIdx] = dStat.getMedian();
		}
		Statistics dStat = new Statistics(sd.locusMedians);
		sd.sampleMedian = dStat.getMedian();
		sd.sampleCoV = dStat.getCoeffOfVariation();
		
		return sd;			
	}

	
	private void writeSampleStats(SampleDepth[] depths) throws AnalysisException {
		String[] locusNames = new String[loci.length];
		for (int i = 0; i < loci.length; i++) {
			locusNames[i] = loci[i].name;
		}
		String[] fieldHeaders = TextUtilities.mergeStringLists(new String[] {"Sample", "SampleMedian", "SampleCoV"}, locusNames);
		TableOutput statOut = null;
		try {
			statOut = new TableOutput(new UncompressedOutputTextStore(outRootFolder, "AllSamples.Medians.tab"), fieldHeaders, 1024 * 1024);
			for (int dIdx = 0; dIdx < depths.length; dIdx++) {
				statOut.newRow();
				statOut.appendValue(depths[dIdx].sample.getName());
				statOut.appendValue(depths[dIdx].sampleMedian);
				statOut.appendValue(depths[dIdx].sampleCoV);
				double[] locusMedians = depths[dIdx].locusMedians;
				for (int lIdx = 0; lIdx < locusMedians.length; lIdx++) {
					statOut.appendValue(locusMedians[lIdx]);
				}
			}
		} finally {
			if (statOut != null) {
				statOut.close();
			}
		}
		try {
			statOut = new TableOutput(new UncompressedOutputTextStore(outRootFolder, "AllSamples.MediansNorm.tab"), fieldHeaders, 1024 * 1024);
			for (int dIdx = 0; dIdx < depths.length; dIdx++) {
				statOut.newRow();
				statOut.appendValue(depths[dIdx].sample.getName());
				statOut.appendValue(depths[dIdx].sampleMedian);
				statOut.appendValue(depths[dIdx].sampleCoV);
				double[] locusMedians = depths[dIdx].locusMedians;
				for (int lIdx = 0; lIdx < locusMedians.length; lIdx++) {
					if (depths[dIdx].sampleMedian == 0.0) {
						statOut.appendBlankValue();
					} else {
						statOut.appendValue(locusMedians[lIdx]/depths[dIdx].sampleMedian);
					}
				}
			}
		} finally {
			if (statOut != null) {
				statOut.close();
			}
		}
	}
	

	private void writeDepths(SampleDepth[] depths, String filename, boolean valuesAsIntegers) throws AnalysisException {
		String[] sampleNames = new String[depths.length];
		for (int i = 0; i < depths.length; i++) {
		    sampleNames[i] = depths[i].sample.getName();
		}
		String[] fieldHeaders = TextUtilities.mergeStringLists(new String[] {"Locus", "Chr", "Pos"}, sampleNames);
		
		TableOutput mergeOut = null;
		try {
			mergeOut = new TableOutput(new GzippedOutputTextStore(outRootFolder, filename), fieldHeaders, 1024 * 1024);
			for (int dpIdx = 0; dpIdx < dpos.length; dpIdx++) {
				DepthPosition dp = dpos[dpIdx];
				mergeOut.newRow();
				mergeOut.appendValue(dp.getName());
				mergeOut.appendValue(dp.getChromosome());
				mergeOut.appendValue(dp.getPos());
				for (int sIdx = 0; sIdx < depths.length; sIdx++) {
					double val = depths[sIdx].depth[dpIdx];
					if (valuesAsIntegers) {
						if (Double.isNaN(val)) {
							mergeOut.appendBlankValue();
						} else {
							mergeOut.appendValue((int)val);
						}
					} else {
						mergeOut.appendValue(val);
					}
				}
			} 
		} finally {
			if (mergeOut != null) {
				mergeOut.close();
			}
		}
	}
	

	private SampleDepth readSampleDepths(Sample sample, DepthPosition[] dpos) throws AnalysisException {
		double[] depth = new double[dpos.length];
		Arrays.fill(depth, Double.NaN);
		
		File sampleFolder = getSampleSubfolder(outRootFolder, sample, true);
		InputTextStore in = null;
		TableInput tif = null;
		try {
			in = new InputTextStore (sampleFolder, sample.getName() + ".depth.tab");
			tif = new TableInput(in);
		} catch (Exception e) {
			log.warn("Skipping sample "+ sample.getName() + ". Could not read depth data: "+e.getMessage());
			return null;
		}
		try {
			String[] inFields;
			int posIdx = 0;
			while ((inFields = tif.getNextValidLine()) != null) {
				DepthPosition currPos = dpos[posIdx];
				String chr = inFields[0];
				int pos = Integer.parseInt(inFields[1]);
				while ((pos != currPos.getPos()) || (!chr.equals(currPos.getChromosome()))) {
					posIdx++;
					if (posIdx < dpos.length) {
						currPos = dpos[posIdx];
					} else {
						throw new AnalysisException("Unknown position in sample file "	+ in.getFile().getAbsolutePath()
						+ " at line " + tif.getLineNumber()	+ ": "+inFields[0]+':'+inFields[1]);
					}
				}
				depth[posIdx] = Integer.parseInt(inFields[2]);
				posIdx++;
			}
		} catch (Exception e) {
			
		} finally {
			tif.close();
		}
		SampleDepth result = new SampleDepth(sample, depth);
		return result;
	}
	
	private DepthLocus[] getLoci (File bedFile) throws AnalysisException {
		HashMap<String,String> chrTable = new HashMap<String,String>();  // Flywheel
		ArrayList<DepthLocus> locusList = new ArrayList<DepthLocus>();
		InputTextStore inStore = new InputTextStore (bedFile);
		TabFileReader r = null;
		try {
			r = new TabFileReader(inStore.getReader());
			String[] inFields;
			while ((inFields = r.getNextValidLine()) != null) {
				// Fields: Chr, startPos, endPos(+1), name
				String chrField = inFields[0];
				String chr = chrTable.get(chrField);
				if (chr == null) {
					chr = new String(chrField);
					chrTable.put(chr,chr);
				}
				int startPos = Integer.parseInt(inFields[1]);
				int endPos = Integer.parseInt(inFields[2]);
				String name = new String(inFields[3]);
				GenomeRegion gr = new GenomeRegion(chr, startPos, endPos);
				locusList.add(new DepthLocus(name, gr));
			}
		} catch (IOException e) {
			throw new AnalysisException("Error processing BED file " + bedFile.getAbsolutePath() +	": "+e);
		} finally {
			try {
				if (r != null) {
					r.close();
				}
			} catch (IOException e) {
				log.info("Error closing BED file: "+e);
			}
		}
		DepthLocus[] loci = locusList.toArray(new DepthLocus[locusList.size()]);
		return loci;
	}
	
	private DepthPosition[] getDepthPositions (DepthLocus[] loci) {
		DepthLocus[] sortedLoci = Arrays.copyOf(loci, loci.length);
		Arrays.sort(sortedLoci);
		int posCount = 0;
		for (DepthLocus locus : sortedLoci) {
			posCount += locus.length;
		}
		DepthPosition[] dpos = new DepthPosition[posCount];
		int dposIdx = 0;
		for (DepthLocus locus : sortedLoci) {
			for (int i = 0 ; i < locus.length; i++) {
				dpos[dposIdx++] = new DepthPosition(locus, i);					
			}
		}
		return dpos;
	}
	
	public class DepthLocus implements Comparable<DepthLocus> {
		String       name;
		GenomeRegion region;
		int          length;
		int          startPosIdx;

		public DepthLocus(String name, GenomeRegion region) {
			this.name = name;
			this.region = region;
			this.length = region.getSize();
		}

		@Override
		public int compareTo(DepthLocus o) {
			return (region.compareTo(o.region));
		}
	}
	
	public class DepthPosition {
		DepthLocus locus;
		int idx;
		
		public DepthPosition(DepthLocus locus, int idx) {
			this.locus = locus;
			this.idx = idx;
		}
		public String getChromosome() {
			return locus.region.getChromosome();
		}
		public int getPos() {
			return (locus.region.getStartPos() + idx);
		}
		public String getName() {
			return locus.name;
		}
	}

	public class SampleDepth {
		Sample   sample;
		double[] depth;
		double   sampleMedian;
		double   sampleCoV;
		double[] locusMedians;
		
		public SampleDepth(Sample sample, double[] depth) {
			this.sample = sample;
			this.depth = depth;
			this.locusMedians = new double[loci.length];
		}
	}
	
	
	public static void main(String[] args) {
		if (args.length < 3) {
			log.error("Usage: org.cggh.bam.depth.DepthAnalysis <sampleListFile> <bedFile> <rootFolder>");
			return;
		}
		File sampleListFile = new File(args[0]);	log.info("SampleListFile: " + sampleListFile.getAbsolutePath());
		File bedFile = new File(args[1]);			log.info("BedFile: " + bedFile.getAbsolutePath());
		File rootFolder = new File(args[2]);		log.info("RootFolder: " + rootFolder.getAbsolutePath());
		
		try {
			Sample[] samples = new SampleList(sampleListFile, false).getSamples();
			DepthAnalysis task = new DepthAnalysis(bedFile, rootFolder);
			task.analyzeAllSamples(samples);
		} catch (Exception e) {
			log.error("Error executing task: " + e);
			return;
		}
		log.info("Exiting");
	}
}
