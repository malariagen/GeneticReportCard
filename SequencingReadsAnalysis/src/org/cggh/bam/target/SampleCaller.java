package org.cggh.bam.target;

import org.cggh.common.counters.LabelCounters;

public class SampleCaller {
	
	private int minStringentCallReadCount;
	private int minLenientCallReadCount;
	
	public SampleCaller(int minStringentCallReadCount, int minLenientCallReadCount) {
		this.minStringentCallReadCount = minStringentCallReadCount;
		this.minLenientCallReadCount = minLenientCallReadCount;
	}
	
	public SampleCall callSample (SampleTargetResult sampleResult) {

		LabelCounters ntAlleleCounters = sampleResult.getNtAlleleCounters();
		Target target = sampleResult.getTarget();
		
		SampleCall ntCall = new SampleCall (ntAlleleCounters, target, minStringentCallReadCount);
		if (ntCall.isMissing()) {
			ntCall = new SampleCall.LenientSampleCall(ntAlleleCounters, target, minLenientCallReadCount);
		}
		return ntCall;
	}
}
