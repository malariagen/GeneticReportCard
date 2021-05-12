package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.counters.*;

public class SampleTargetResult {

	private Target          target;
	private Sample          sample;
	private SampleCall      ntCall;
	private AminoSampleCall aaCall;
	private LabelCounters   ntAlleleCounters;
	private int             lowQualityCount;

	public SampleTargetResult(Target target, Sample sample, SampleCall ntCall, AminoSampleCall aaCall, LabelCounters ntCounters, int lowQualityCount) {
		this.target = target;
		this.sample = sample;
		this.ntCall = ntCall;
		this.aaCall = aaCall;
		this.ntAlleleCounters = ntCounters;
		this.lowQualityCount = lowQualityCount;
	}

	public SampleTargetResult(Target target, Sample sample) {
		this(target, sample, SampleCall.makeMissingCall(), AminoSampleCall.makeMissingCall(), new LabelCounters(), 0);
	}

	public Target getTarget() {
		return target;
	}

	public Sample getSample() {
		return sample;
	}

	public SampleCall getNtCall() {
		return ntCall;
	}

	public AminoSampleCall getAminoCall() {
		return aaCall;
	}

	public int getLowQualityCount() {
		return lowQualityCount;
	}
	
	public LabelCounters getNtAlleleCounters() {
		return ntAlleleCounters;
	}
	
	public void setNtAlleleCounters(LabelCounters ntAlleleCounters) {
		this.ntAlleleCounters = ntAlleleCounters;
	}
	
	public void filterAlleleCounters (LabelCounterFilter filter) {
		ntAlleleCounters.filterCounters(filter);
	}
	
	public void cleanupTargetAlleles () {
		// Clean up the target sequence
		filterAlleleCounters (new LabelCounterFilterByLabelSubstring("N"));  // Remove sequences with undetermined nucleotide
		filterAlleleCounters (new LabelCounterFilterByMinCount(2));          // Remove Singletons			
	}
}
