package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import java.util.*;


public class TargetLocus extends Locus {

	private Target[] targets;
	
	public TargetLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, boolean analyzeUnmappedReads, Target[] targets) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, analyzeUnmappedReads);
		this.targets = targets;
	}

	public static TargetLocus[] parseLocusConfig(Properties configProperties, String propPrefix) throws AnalysisException {
		return new TargetLocusConfig(configProperties).parseLoci(propPrefix);
	}

	public Target[] getTargets() {
		return targets;
	}

	public static class TargetLocusConfig extends LocusConfig {

		public TargetLocusConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public TargetLocus[] parseLoci (String propPrefix) throws AnalysisException {
			
			Locus[] baseLoci = super.parseLoci(propPrefix);
			
			TargetLocus[] loci = new TargetLocus[baseLoci.length];
			for (int lIdx = 0; lIdx < loci.length; lIdx++) {
				String propLocusPrefix = propPrefix+"locus."+baseLoci[lIdx].getName();
				String locusChromosome = baseLoci[lIdx].getChromosome();
				Target[] targets = Target.parseTargetConfig(configProperties, propLocusPrefix, locusChromosome);
				Locus baseLocus = baseLoci[lIdx];
				loci[lIdx] = new TargetLocus(baseLocus.getName(), baseLocus.getReadSearchIntervalCoords(), baseLocus.getAnchors(), baseLocus.getAnalyzeUnmappedReads(), targets);
			}
			return loci;
		}
	}
}
