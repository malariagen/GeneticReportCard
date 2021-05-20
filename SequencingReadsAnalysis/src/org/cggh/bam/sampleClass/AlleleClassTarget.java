package org.cggh.bam.sampleClass;

import org.cggh.bam.target.Target;
import org.cggh.common.exceptions.*;
import java.util.*;


public class AlleleClassTarget extends Target {
	
	private ClassAllele[] alleles;
	
	public AlleleClassTarget (String name, String[] targetCoords, boolean isReverse, ClassAllele[] alleles) throws AnalysisException {
		super (name, targetCoords, isReverse);
		this.alleles = alleles;
	}

	public ClassAllele[] getAlleles() {
		return alleles;
	}

	public static class ClassAllele {
		private String name;
		private String[] sequences;
		
		public ClassAllele(String name, String[] sequences) {
			this.name = name;
			this.sequences = sequences;
		}

		public String getName() {
			return name;
		}

		public String[] getSequences() {
			return sequences;
		}
	}
	
	
	public static AlleleClassTarget[] parseTargetConfig(Properties configProperties, String propLocusPrefix, String locusChromosome) throws AnalysisException {
		return new AlleleClassTargetConfig(configProperties).parseTargets(propLocusPrefix, locusChromosome);
	}
	
	public static class AlleleClassTargetConfig extends TargetConfig {

		public AlleleClassTargetConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public AlleleClassTarget[] parseTargets (String propLocusPrefix, String locusChromosome) throws AnalysisException {
						
			Target[] baseTargets = super.parseTargets(propLocusPrefix, locusChromosome);
			
			AlleleClassTarget[] targets = new AlleleClassTarget[baseTargets.length];
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				String targetPrefix = propLocusPrefix+".target."+baseTargets[tIdx].getName();
				String[] allelesStr = getStringListProperty(targetPrefix+".alleles");
				ClassAllele[] alleles = new ClassAllele[allelesStr.length];
				for (int aIdx = 0; aIdx < alleles.length; aIdx++) {
					String[] parts = allelesStr[aIdx].split("@");
					String alleleName = parts[0].trim();
					String[] alleleSequences = parts[1].trim().split("\\|");
					alleles[aIdx] = new ClassAllele (alleleName, alleleSequences);
				}
				Target baseTarget = baseTargets[tIdx];
				targets[tIdx] = new AlleleClassTarget(baseTarget.getName(), baseTarget.getTargetCoords(), baseTarget.isReverse(), alleles);
			}
			return targets;
		}
	}
}


