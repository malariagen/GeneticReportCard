package org.cggh.bam.target.amino;

import org.cggh.common.counters.*;
import org.cggh.common.sequence.*;

public class AminoAlleleCounter extends LabelCounter {
	
	public String ntSequence;
	public String aaSequence;
	public String fullLabel;
	
	public AminoAlleleCounter (String ntSequence) {
		super (ntSequence);
		this.ntSequence = ntSequence;
		this.aaSequence = SequenceUtilities.translateNtSequence(ntSequence);
		this.fullLabel = ntSequence+"["+aaSequence+"]";
	}
	
	public String getSummary() {
		return appendCountToLabel(fullLabel);
	}	
}

