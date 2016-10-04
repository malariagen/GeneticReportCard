package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;

public class Target {


	String       name;
	GenomeRegion targetRegion;
	boolean      isReverse;
	String       targetNtSeq;
	String       targetAaSeq;
	
	public Target (String name, String targetCoords, boolean isReverse) throws AnalysisException {
		this.name = name;
		this.targetRegion = GenomeRegion.parseRegion(targetCoords);
		this.isReverse = isReverse;

		Sequence chrSeq = ReferenceGenome.getChrSequence(targetRegion.getChromosome());
		int startIdx = targetRegion.getStartPos() - 1;
		int endIdx = targetRegion.getStopPos();
	    targetNtSeq = chrSeq.getData().substring(startIdx, endIdx);
	    if (isReverse) {
		    targetNtSeq = SequenceUtilities.getReverseComplementSequence(targetNtSeq);
	    }
	    targetAaSeq = SequenceUtilities.translateNtSequence(targetNtSeq);
	}
	
	public String getName() {
		return name;
	}

	public GenomeRegion getTargetRegion() {
		return targetRegion;
	}

	public boolean isReverse() {
		return isReverse;
	}

	public String getTargetNtSeq() {
		return targetNtSeq;
	}

	public String getTargetAaSeq() {
		return targetAaSeq;
	}
}


