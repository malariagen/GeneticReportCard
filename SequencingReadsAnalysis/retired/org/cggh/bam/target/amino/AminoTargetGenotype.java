package org.cggh.bam.target.amino;

import org.cggh.bam.target.*;
import org.cggh.common.sequence.*;

public class AminoTargetGenotype extends TargetGenotype {
	
	protected String      aminoGenotype;
	
	public AminoTargetGenotype(String ntGenotype, AminoTarget target) {
		super(ntGenotype, target);
		aminoGenotype = (ntGenotype == null) ? null : SequenceUtilities.translateNtSequence(ntGenotype);
	}
	
	public String getAminoGenotype () {
		return aminoGenotype;
	}
	
	public boolean isReferenceAllele () {
	    if ((aminoGenotype == null) || (target == null)) {
	    	return false;
	    }
	    return (aminoGenotype.equals(((AminoTarget)target).getAminoTargetRefSeq()));
	}
}

