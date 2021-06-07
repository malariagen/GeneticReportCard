package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;


public class TargetLocus extends Locus {

	protected Target[] targets;
	
	public TargetLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, analyzeUnmappedReads);
	}

	public void setTargets(Target[] targets) {
		this.targets = targets;
		for (int i = 0; i < targets.length; i++) {
			targets[i].setLocus(this);
		}
	}

	public Target[] getTargets() {
		return targets;
	}
}
