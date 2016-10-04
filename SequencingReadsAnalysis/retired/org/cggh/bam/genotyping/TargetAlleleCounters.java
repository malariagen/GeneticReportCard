package org.cggh.bam.genotyping;

import org.cggh.common.counters.LabelCounter;
import org.cggh.common.counters.LabelCounters;

public class TargetAlleleCounters extends LabelCounters {
	public LabelCounter createCounter (String ntSequence) {
		return new TargetAlleleCounter(ntSequence);
	}
}
