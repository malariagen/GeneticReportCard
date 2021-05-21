package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;

public class Target {

	private String         name;
	private String[]       targetCoords;
	private GenomeRegion[] targetRegions;
	private boolean        isReverse;
	private String         targetRefSeq;

	
	public Target (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		this.name = name;
		this.targetCoords = targetCoords;
		targetRegions = new GenomeRegion[targetCoords.length];
		for (int i = 0; i < targetCoords.length; i++) {
			targetRegions[i] = GenomeRegion.parseRegion(targetCoords[i]);
		}
		this.isReverse = isReverse;

		Sequence chrSeq = ReferenceGenome.getChrSequence(targetRegions[0].getChromosome());
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < targetRegions.length; i++) {
			GenomeRegion targetRegion = targetRegions[i];
			int startIdx = targetRegion.getStartPos() - 1;
			int endIdx = targetRegion.getStopPos();
			sb.append(chrSeq.getData().substring(startIdx, endIdx));
		}
		targetRefSeq = sb.toString();
	    if (isReverse) {
	    	targetRefSeq = SequenceUtilities.getReverseComplementSequence(targetRefSeq);
	    }
	}
	
	public String getName() {
		return name;
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

	public String getTargetRefSeq() {
		return targetRefSeq;
	}
	
}


