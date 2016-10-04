package org.cggh.bam.target.amino;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import java.util.*;

public class AminoTargetGenotyper extends TargetGenotyper {
	
	protected AminoTarget target;
	
	public AminoTargetGenotyper (AminoTarget target) {
		super(target);
		this.target = target;
	}
		
	@Override
	public TargetGenotype[] extractTargetGenotypes(MappedRead[] reads) {
		HashMap<String,AminoTargetGenotype> genoTable = new HashMap<String,AminoTargetGenotype>(reads.length);
		TargetGenotype[] baseGenos = super.extractTargetGenotypes(reads);
		TargetGenotype[] targetGenos = new TargetGenotype[reads.length];
		for (int rIdx = 0; rIdx < reads.length; rIdx++) {
			String ntSequence = baseGenos[rIdx].getNtGenotype();
			if (ntSequence == null) {
				targetGenos[rIdx] = baseGenos[rIdx];
				continue;
			}
			AminoTargetGenotype aaGeno = genoTable.get(ntSequence);
			if (aaGeno == null) {
				aaGeno = new AminoTargetGenotype (ntSequence, target);
				genoTable.put(ntSequence, aaGeno);
			}
			targetGenos[rIdx] = aaGeno;
		}
		return targetGenos;	
	}
}
