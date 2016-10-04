package org.cggh.bam.target.amino;

import org.cggh.bam.target.Target;
import org.cggh.common.exceptions.*;
import org.cggh.common.sequence.SequenceUtilities;

import java.util.*;


public class AminoTarget extends Target {

	private String aminoTargetRefSeq;
	
	public AminoTarget (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		super (name, targetCoords, isReverse);
		aminoTargetRefSeq = SequenceUtilities.translateNtSequence(getTargetRefSeq());
	}
	
	public String getAminoTargetRefSeq() {
		return aminoTargetRefSeq;
	}
	
	
	public static AminoTarget[] parseTargetConfig(Properties configProperties, String propLocusPrefix, String locusChromosome) throws AnalysisException {
		return new AminoTargetConfig(configProperties).parseTargets(propLocusPrefix, locusChromosome);
	}
	
	public static class AminoTargetConfig extends TargetConfig {

		public AminoTargetConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public AminoTarget[] parseTargets (String propLocusPrefix, String locusChromosome) throws AnalysisException {		
			Target[] baseTargets = super.parseTargets(propLocusPrefix, locusChromosome);
			AminoTarget[] targets = new AminoTarget[baseTargets.length];
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				Target baseTarget = baseTargets[tIdx];
				targets[tIdx] = new AminoTarget(baseTarget.getName(), baseTarget.getTargetCoords(), baseTarget.isReverse());
			}
			return targets;
		}
	}
}


