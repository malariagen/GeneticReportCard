package org.cggh.bam.sampleClass;

import org.cggh.bam.MappedRead;
import org.cggh.bam.sampleClass.AlleleClassGenotype.*;
import org.cggh.bam.sampleClass.AlleleClassTarget.*;
import org.cggh.bam.target.*;

public class AlleleClassGenotyper extends TargetGenotyper {
	
	protected AlleleClassTarget           target;
	protected ClassAllele[]               alleles;
	protected AlleleClassGenotype[] definedGenos;
	
	public AlleleClassGenotyper (AlleleClassTarget target) {
		super(target);
		this.target = target;
		alleles = target.getAlleles();
		definedGenos = new AlleleClassGenotype[alleles.length];
		for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
			definedGenos[aIdx] = new AlleleClassGenotype((AlleleClassTarget)target, aIdx);
		}
	}
	
	public TargetGenotype[] extractTargetNtGenotypes (MappedRead[] reads) {
		TargetGenotype[] baseGenos = super.extractTargetGenotypes(reads);
		TargetGenotype[] targetGenos = new TargetGenotype[reads.length];
		for (int rIdx = 0; rIdx < reads.length; rIdx++) {
			String geno = baseGenos[rIdx].getNtGenotype();
			if (geno == null) {
				targetGenos[rIdx] = baseGenos[rIdx];
				continue;
			}
			for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
				String[] sequences = alleles[aIdx].getSequences();
				for (int sIdx = 0; sIdx < sequences.length; sIdx++) {
					if (geno.equals(sequences[sIdx])) {
						targetGenos[rIdx] = definedGenos[aIdx];
						break;
					}
				}
				if (targetGenos[rIdx] != null) {
					 break;
				}
			}
			if (targetGenos[rIdx] == null) {
				targetGenos[rIdx] = new UnlistedAlleleTargetGenotype(geno, target);
			}
		}
		return targetGenos;	
	}
}
