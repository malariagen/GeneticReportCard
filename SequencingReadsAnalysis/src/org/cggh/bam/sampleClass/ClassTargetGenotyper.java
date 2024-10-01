package org.cggh.bam.sampleClass;

import org.cggh.bam.BaseAnalysisConfig;
import org.cggh.bam.Read;
import org.cggh.bam.sampleClass.ClassTarget.*;
import org.cggh.bam.target.*;
import org.cggh.common.util.TextUtilities;

public class ClassTargetGenotyper extends TargetGenotyper {
	
	protected ClassTarget     target;
	protected ClassAllele[]         alleles;
	protected ClassTargetGenotype[] definedGenos;
	
	public ClassTargetGenotyper (ClassTarget target, BaseAnalysisConfig config) {
		super(target, config);
		this.target = target;
		alleles = target.getAlleles();
		definedGenos = new ClassTargetGenotype[alleles.length];
		for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
			definedGenos[aIdx] = new ClassTargetGenotype((ClassTarget)target, aIdx);
		}
	}
	
	public TargetGenotype[] extractTargetNtGenotypes (Read[] reads) {
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
	
	public static class ClassTargetGenotype extends TargetGenotype {
		ClassTarget target;
		int alleleIdx;
		
		public ClassTargetGenotype(ClassTarget target, int alleleIdx) {
			super(TextUtilities.stringArrayToString (target.getAlleles()[alleleIdx].getSequences()), target);
			this.target = target;
			this.alleleIdx = alleleIdx;
		}
		public ClassAllele getAllele () {
			return target.getAlleles()[alleleIdx];
		}
	}
	
	public static class UnlistedAlleleTargetGenotype extends TargetGenotype {
		public UnlistedAlleleTargetGenotype (String seq, ClassTarget target) {
			super(seq, target);
		}
	};
}
