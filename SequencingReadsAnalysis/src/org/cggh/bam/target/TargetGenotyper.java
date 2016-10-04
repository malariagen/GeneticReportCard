package org.cggh.bam.target;

import java.util.HashMap;

import org.cggh.bam.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;

public class TargetGenotyper {
	
	protected Target target;
	
	public TargetGenotyper (Target target) {
		this.target = target;
	}
	
	public TargetGenotype[] extractTargetGenotypes (MappedRead[] reads) {
		GenomeRegion[] tRegions = target.getTargetRegions();
		HashMap<String,TargetGenotype> genoTable = new HashMap<String,TargetGenotype>(reads.length);
		TargetGenotype[] targetGenos = new TargetGenotype[reads.length];
		StringBuffer sb = new StringBuffer();
		
		nextRead:
		for (int rIdx = 0; rIdx < reads.length; rIdx++) {
			MappedRead r = reads[rIdx];
			int rStartPos = r.getStartPos();
			int rEndPos = rStartPos + r.getSequence().length() - 1;
			
			sb.setLength(0);
			for (int trIdx = 0; trIdx < tRegions.length; trIdx++) {
				GenomeRegion tRegion = tRegions[trIdx];
				int tStartPos = tRegion.getStartPos();
				int tEndPos = tRegion.getStopPos();
				int tLen = 1 + tEndPos - tStartPos;

				if ((rStartPos > tStartPos) || (rEndPos < tEndPos)) {
					targetGenos[rIdx] = new TargetGenotype.NoCoverageTargetGenotype();
					continue nextRead;
				}
				
				int tStartOffset = tStartPos - rStartPos;
				String trSeq = r.getSequence().substring(tStartOffset, tStartOffset+tLen);
				if (trSeq.contains("N")) {
					targetGenos[rIdx] = new TargetGenotype.NoCoverageTargetGenotype();
					continue nextRead;
				}
				
				int minQ = 1000;
				for (int j = 0; j < tLen; j++) {
					int q = r.getPhredScore(tStartOffset+j);
					if (q < minQ) {
						minQ = q;
					}
				}
				if (minQ < MappedRead.MIN_PHRED_SCORE) {
					targetGenos[rIdx] = new TargetGenotype.LowQualityTargetGenotype();
					continue nextRead;
				}
				sb.append(trSeq);
			}
			String ntSequence = sb.toString();
			
			// If the gene is negative-strand, reverse the sequence
			if (target.isReverse()) {
				ntSequence = SequenceUtilities.getReverseComplementSequence(ntSequence);					
			}
			
			// Get instance from table (flywheel)
			TargetGenotype ntGeno = genoTable.get(ntSequence);
			if (ntGeno == null) {
				ntGeno = new TargetGenotype (ntSequence, target);
				genoTable.put(ntSequence, ntGeno);
			}
			targetGenos[rIdx] = ntGeno;
		}
		return targetGenos;	
	}
}



