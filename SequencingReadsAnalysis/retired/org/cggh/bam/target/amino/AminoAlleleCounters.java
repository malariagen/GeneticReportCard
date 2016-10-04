package org.cggh.bam.target.amino;

import org.cggh.common.counters.LabelCounter;
import org.cggh.common.counters.LabelCounters;

public class AminoAlleleCounters extends LabelCounters {
	public LabelCounter createCounter (String ntSequence) {
		return new AminoAlleleCounter(ntSequence);
	}
}
