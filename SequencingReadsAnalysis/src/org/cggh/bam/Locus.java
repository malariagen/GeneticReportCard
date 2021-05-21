package org.cggh.bam;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.genome.GenomeRegion;

public class Locus {

	String       name;	
	String       readSearchIntervalCoords;
	GenomeRegion readSearchInterval;
	Anchor[]     anchors;
	boolean      analyzeUnmappedReads;

	public Locus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads) throws AnalysisException {
		this.name = name;
		this.readSearchIntervalCoords = readSearchIntervalCoords;
		this.readSearchInterval = GenomeRegion.parseRegion(readSearchIntervalCoords);
		this.anchors = anchors;
		this.analyzeUnmappedReads = analyzeUnmappedReads;
	}

	public String getName() {
		return name;
	}

	public String getChromosome() {
		return readSearchInterval.getChromosome();
	}

	public String getReadSearchIntervalCoords() {
		return readSearchIntervalCoords;
	}

	public GenomeRegion getReadSearchInterval() {
		return readSearchInterval;
	}

	public Anchor[] getAnchors() {
		return anchors;
	}

	public boolean getAnalyzeUnmappedReads() {
		return analyzeUnmappedReads;
	}
}
