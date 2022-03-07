package org.cggh.bam.sampleClass;

import org.cggh.common.counters.*;

public class SampleTargetResult {

	private ClassTarget    target;
	private int            totalReadCount;
	private LabelCounter[] alleleCounters;
	private LabelCounter[] unlistedAlleleCounters;
	private String[]       classTargetCalls;
	private String         call;

	public SampleTargetResult(ClassTarget target, String[] classTargetCalls, LabelCounter[] alleleCounters, LabelCounter[] unlistedAlleleCounters) {
		this.target = target;
		this.classTargetCalls = classTargetCalls;
		this.alleleCounters = alleleCounters;
		this.unlistedAlleleCounters = unlistedAlleleCounters;
		totalReadCount = 0;
		for (int aIdx = 0; aIdx < alleleCounters.length; aIdx++) {
			totalReadCount += alleleCounters[aIdx].getCount();
		}
		for (int aIdx = 0; aIdx < unlistedAlleleCounters.length; aIdx++) {
			totalReadCount += unlistedAlleleCounters[aIdx].getCount();
		}
		for (int tcIdx = 0; tcIdx < classTargetCalls.length; tcIdx++) {
			String alleleLabel = classTargetCalls[tcIdx];
			if ("-".equals(alleleLabel)) {
				continue;
			}
			call = (call == null) ? alleleLabel : call+","+alleleLabel;
		}
	}
	
	public ClassTarget getTarget() {
		return target;
	}

	public String[] getClassTargetCalls() {
		return classTargetCalls;
	}

	public String getCall() {
		return call;
	}

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public LabelCounter[] getAlleleCounters() {
		return alleleCounters;
	}

	public LabelCounter[] getUnlistedAlleleCounters() {
		return unlistedAlleleCounters;
	}
}
