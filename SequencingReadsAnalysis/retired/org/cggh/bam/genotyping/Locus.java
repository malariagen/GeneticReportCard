package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.genome.GenomeRegion;

public class Locus {

	String   name;
	//String   readSearchIntervalCoords;	
	GenomeRegion readSearchInterval;
	Target[] targets;
	Anchor[] anchors;

	public Locus(String name, String readSearchIntervalCoords, Target[] targets, Anchor[] anchors) throws AnalysisException {
		this.name = name;
		this.readSearchInterval = GenomeRegion.parseRegion(readSearchIntervalCoords);
		this.targets = targets;
		this.anchors = anchors;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public GenomeRegion getReadSearchInterval() {
		return readSearchInterval;
	}

	public void setReadSearchInterval(GenomeRegion readSearchInterval) {
		this.readSearchInterval = readSearchInterval;
	}

	public Target[] getTargets() {
		return targets;
	}

	public void setTargets(Target[] targets) {
		this.targets = targets;
	}

	public Anchor[] getAnchors() {
		return anchors;
	}

	public void setAnchors(Anchor[] anchors) {
		this.anchors = anchors;
	}		
}
