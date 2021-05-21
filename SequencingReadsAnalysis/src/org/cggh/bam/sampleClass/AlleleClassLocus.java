package org.cggh.bam.sampleClass;

import org.cggh.bam.Anchor;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.*;


public class AlleleClassLocus extends TargetLocus {
	
	public AlleleClassLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads, AlleleClassTarget[] targets) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, analyzeUnmappedReads, targets);
	}

	public AlleleClassTarget[] getTargets() {
		return (AlleleClassTarget[])targets;
	}
}
