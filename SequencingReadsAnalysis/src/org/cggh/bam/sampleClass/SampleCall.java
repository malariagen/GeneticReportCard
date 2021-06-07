package org.cggh.bam.sampleClass;

import org.cggh.bam.*;

public class SampleCall {
	
	private Sample               sample;
	private String               classCall;
	private SampleTargetResult[] targetResults;

	public SampleCall(Sample sample, String classCall, SampleTargetResult[] targetResults) {
		this.sample = sample;
		this.classCall = classCall;
		this.targetResults = targetResults;
	}
	
	public Sample getSample() {
		return sample;
	}

	public String getClassCall() {
		return classCall;
	}

	public SampleTargetResult[] getTargetResults() {
		return targetResults;
	}
}
