package org.cggh.bam.sampleClass;

import org.cggh.bam.Anchor;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.*;


public class ClassLocus extends TargetLocus {
	
	public ClassLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, analyzeUnmappedReads);
	}

	public void setTargets(ClassTarget[] targets) {
		super.setTargets(targets);
	}

	public ClassTarget[] getTargets() {
		return (ClassTarget[])targets;
	}
}
