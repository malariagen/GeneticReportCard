package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;

public class AlignmentTarget extends Target{

	private String         targetRefSeq;

	public AlignmentTarget (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		super (name, targetCoords, isReverse);

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
	
	public String getTargetRefSeq() {
		return targetRefSeq;
	}
}


