package org.cggh.bam;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.genome.GenomeRegion;

public class Locus {

	String         name;	
	String         readSearchIntervalCoords;
	GenomeRegion[] readSearchIntervals;
	Anchor[]       anchors;
	boolean        analyzeUnmappedReads;

	public Locus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads) throws AnalysisException {
		this.name = name;
		
		this.readSearchIntervalCoords = readSearchIntervalCoords;
		String[] intervalStrs = readSearchIntervalCoords.split(",");
		this.readSearchIntervals = new GenomeRegion[intervalStrs.length];
		for (int idx = 0; idx < intervalStrs.length; idx++) {
			this.readSearchIntervals[idx] = GenomeRegion.parseRegion(intervalStrs[idx]);
		}
		
		this.anchors = anchors;
		this.analyzeUnmappedReads = analyzeUnmappedReads;
	}

	public String getName() {
		return name;
	}

	public String getReadSearchIntervalCoords() {
		return readSearchIntervalCoords;
	}

	public GenomeRegion[] getReadSearchIntervals() {
		return readSearchIntervals;
	}

	// Convenience function for alignment-based tasks
	public GenomeRegion getReadSearchInterval() {
		return readSearchIntervals[0];
	}

	public Anchor[] getAnchors() {
		return anchors;
	}

	public boolean getAnalyzeUnmappedReads() {
		return analyzeUnmappedReads;
	}
}
