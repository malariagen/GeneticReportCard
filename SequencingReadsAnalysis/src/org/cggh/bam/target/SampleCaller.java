package org.cggh.bam.target;

import org.cggh.common.counters.LabelCounters;

public class SampleCaller {
	
	private int minStringentCallReadCount;
	private int minStringentAlleleReadCount;
	private int minLenientCallReadCount;
	private int minLenientAlleleReadCount;
	
	public SampleCaller(int minStringentCallReadCount, int minStringentAlleleReadCount, int minLenientCallReadCount, int minLenientAlleleReadCount) {
		this.minStringentCallReadCount = minStringentCallReadCount;
		this.minStringentAlleleReadCount = minStringentAlleleReadCount;
		this.minLenientCallReadCount = minLenientCallReadCount;
		this.minLenientAlleleReadCount = minLenientAlleleReadCount;
	}
	
	public SampleCall callSample (SampleTargetResult sampleResult) {
		LabelCounters ntAlleleCounters = sampleResult.getNtAlleleCounters();
		Target target = sampleResult.getTarget();
		
		SampleCall ntCall = new SampleCall (ntAlleleCounters, target, minStringentCallReadCount, minStringentAlleleReadCount);
		if (ntCall.isMissing()) {
			ntCall = new SampleCall.LenientSampleCall(ntAlleleCounters, target, minLenientCallReadCount, minLenientAlleleReadCount);
		}
		return ntCall;
	}
}
