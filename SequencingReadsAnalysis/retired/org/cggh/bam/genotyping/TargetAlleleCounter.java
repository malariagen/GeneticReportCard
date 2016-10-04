package org.cggh.bam.genotyping;

import org.cggh.common.counters.*;
import org.cggh.common.sequence.*;

public class TargetAlleleCounter extends LabelCounter {
	
	public String ntSequence;
	public String aaSequence;
	public String fullLabel;
	
	public TargetAlleleCounter (String ntSequence) {
		super (ntSequence);
		this.ntSequence = ntSequence;
		this.aaSequence = SequenceUtilities.translateNtSequence(ntSequence);
		this.fullLabel = ntSequence+"["+aaSequence+"]";
	}
	
	public String getSummary() {
		return appendCountToLabel(fullLabel);
	}	
}

