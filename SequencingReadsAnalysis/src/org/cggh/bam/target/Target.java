package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;

public class Target {

	protected String         name;
	protected TargetLocus    locus;
	protected String[]       targetCoords;
	protected GenomeRegion[] targetRegions;
	protected boolean        isReverse;

	public Target (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		this.name = name;
		this.targetCoords = targetCoords;
		targetRegions = new GenomeRegion[targetCoords.length];
		for (int i = 0; i < targetCoords.length; i++) {
			targetRegions[i] = GenomeRegion.parseRegion(targetCoords[i]);
		}
		this.isReverse = isReverse;
	}
	
	public void setLocus(TargetLocus locus) {
		this.locus = locus;
	}

	public String getName() {
		return name;
	}

	public Locus getLocus() {
		return locus;
	}

	public String[] getTargetCoords() {
		return targetCoords;
	}

	public GenomeRegion[] getTargetRegions() {
		return targetRegions;
	}

	public boolean isReverse() {
		return isReverse;
	}
}


