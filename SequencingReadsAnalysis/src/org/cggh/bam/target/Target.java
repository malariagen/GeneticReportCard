package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import org.cggh.common.sequence.*;
import java.util.*;

public class Target {

	private String         name;
	private String[]       targetCoords;
	private GenomeRegion[] targetRegions;
	private boolean        isReverse;
	private String         targetRefSeq;

	
	public Target (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		this.name = name;
		this.targetCoords = targetCoords;
		targetRegions = new GenomeRegion[targetCoords.length];
		for (int i = 0; i < targetCoords.length; i++) {
			targetRegions[i] = GenomeRegion.parseRegion(targetCoords[i]);
		}
		this.isReverse = isReverse;

		Sequence chrSeq = ReferenceGenome.getChrSequence(targetRegions[0].getChromosome());
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < targetRegions.length; i++) {
			GenomeRegion targetRegion = targetRegions[i];
			int startIdx = targetRegion.getStartPos() - 1;
			int endIdx = targetRegion.getStopPos();
			sb.append(chrSeq.getData().substring(startIdx, endIdx));
		}
		targetRefSeq = sb.toString();
	    if (isReverse) {
	    	targetRefSeq = SequenceUtilities.getReverseComplementSequence(targetRefSeq);
	    }
	}
	
	public String getName() {
		return name;
	}

	public String[] getTargetCoords() {
		return targetCoords;
	}

	public GenomeRegion[] getTargetRegions() {
		return targetRegions;
	}

	public boolean isReverse() {
		return isReverse;
	}

	public String getTargetRefSeq() {
		return targetRefSeq;
	}
	
	public static Target[] parseTargetConfig(Properties configProperties, String propLocusPrefix, String locusChromosome) throws AnalysisException {
		return new TargetConfig(configProperties).parseTargets(propLocusPrefix, locusChromosome);
	}
	
	public static class TargetConfig extends BaseConfig {

		public TargetConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public Target[] parseTargets (String propLocusPrefix, String locusChromosome) throws AnalysisException {

			String[] targetsStr = getStringListProperty(propLocusPrefix+".targets");
			Target[] targets = new Target[targetsStr.length];
			for (int idx = 0; idx < targets.length; idx++) {
				String[] parts = targetsStr[idx].split("@");
				String coordStr = parts[1].trim();
				boolean isReverse = false;
				if (coordStr.startsWith("-")) {
					isReverse = true;
					coordStr = coordStr.substring(1);  // Remove -
				}
				String[] coordStrs = coordStr.split("&");
				for (int j = 0; j < coordStrs.length; j++) {
					coordStrs[j] = locusChromosome+":"+coordStrs[j];
				}
				targets[idx] = new Target(parts[0], coordStrs, isReverse);
			}
			return targets;
		}
	}
}


