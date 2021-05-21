package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;


public class TargetLocus extends Locus {

	protected Target[] targets;
	
	public TargetLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads, Target[] targets) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, analyzeUnmappedReads);
		this.targets = targets;
	}

	public Target[] getTargets() {
		return targets;
	}
}
