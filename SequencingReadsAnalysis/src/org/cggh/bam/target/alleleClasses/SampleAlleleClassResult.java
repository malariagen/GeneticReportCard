package org.cggh.bam.target.alleleClasses;

import org.cggh.common.counters.*;

public class SampleAlleleClassResult {

	private AlleleClassTarget target;
	private int               totalReadCount;
	private LabelCounter[]    alleleSetCounters;
	private LabelCounter[]    unlistedAlleleCounters;

	public SampleAlleleClassResult(AlleleClassTarget target, LabelCounter[] alleleSetCounters, LabelCounter[] unlistedAlleleCounters) {
		this.target = target;
		this.alleleSetCounters = alleleSetCounters;
		this.unlistedAlleleCounters = unlistedAlleleCounters;
		totalReadCount = 0;
		for (int aIdx = 0; aIdx < alleleSetCounters.length; aIdx++) {
			totalReadCount += alleleSetCounters[aIdx].getCount();
		}
		for (int aIdx = 0; aIdx < unlistedAlleleCounters.length; aIdx++) {
			totalReadCount += unlistedAlleleCounters[aIdx].getCount();
		}
	}

	public AlleleClassTarget getTarget() {
		return target;
	}

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public LabelCounter[] getAlleleSetCounters() {
		return alleleSetCounters;
	}

	public LabelCounter[] getUnlistedAlleleCounters() {
		return unlistedAlleleCounters;
	}
}
