package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.counters.*;

public class SampleTargetResult {

	private Target        target;
	private Sample        sample;
	private LabelCounters ntAlleleCounters;
	

	public SampleTargetResult(Target target, Sample sample) {
		this.target = target;
		this.sample = sample;
		this.ntAlleleCounters = new LabelCounters();
	}

	public SampleTargetResult(Target target, Sample sample, LabelCounters ntCounters) {
		this.target = target;
		this.sample = sample;
		this.ntAlleleCounters = ntCounters;
	}

	public Target getTarget() {
		return target;
	}

	public Sample getSample() {
		return sample;
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
