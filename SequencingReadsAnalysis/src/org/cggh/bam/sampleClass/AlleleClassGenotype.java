package org.cggh.bam.sampleClass;

import org.cggh.bam.sampleClass.AlleleClassTarget.*;
import org.cggh.bam.target.TargetGenotype;
import org.cggh.common.util.*;

public class AlleleClassGenotype extends TargetGenotype {
	
	AlleleClassTarget target;
	int alleleIdx;
	
	public AlleleClassGenotype(AlleleClassTarget target, int alleleIdx) {
		super(TextUtilities.stringArrayToString (target.getAlleles()[alleleIdx].getSequences()), target);
		this.target = target;
		this.alleleIdx = alleleIdx;
	}
	
	public ClassAllele getAllele () {
		return target.getAlleles()[alleleIdx];
	}

	public static class UnlistedAlleleTargetGenotype extends TargetGenotype {
		public UnlistedAlleleTargetGenotype (String seq, AlleleClassTarget target) {
			super(seq, target);
		}
	};
}

