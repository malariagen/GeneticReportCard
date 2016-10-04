package org.cggh.bam.genotyping;

import org.cggh.common.counters.*;

public class AlignmentTarget {
	
	Target               target;
	int                  alleleCount = 0;
	TargetAlleleCounters alleleCounters    = new TargetAlleleCounters();
	TargetAlleleCounters rawAlleleCounters = new TargetAlleleCounters();
	SampleCall           stringentCall;
	SampleCall           lenientCall;
	
	public AlignmentTarget(Target target) {
		this.target = target;
	}

	public void filterAlleleCounters (LabelCounterFilter filter) {
		alleleCounters.filterCounters(filter);
	}
	
	public void cleanupTargetAlleles () {
		// Clean up the target sequence
		filterAlleleCounters (new LabelCounterFilterByLabelSubstring("N"));  // Remove sequences with undetermined nucleotide
		filterAlleleCounters (new LabelCounterFilterByMinCount(2));          // Remove Singletons			
	}
	
	public void callGenotypes () {
		stringentCall = new SampleCall (alleleCounters, target, 5, 2);
		lenientCall = new SampleCall (alleleCounters, target, 2, 2);
	}
}
