package org.cggh.bam.target.alleleClasses;

import org.cggh.bam.Anchor;
import org.cggh.bam.Locus;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.*;
import java.util.*;


public class AlleleClassLocus extends TargetLocus {

	private AlleleClassTarget[] targets;
	
	public AlleleClassLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, AlleleClassTarget[] targets) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors, targets);
		this.targets = targets;
	}

	public static AlleleClassLocus[] parseLocusConfig(Properties configProperties, String propPrefix) throws AnalysisException {
		return new AlleleClassLocusConfig(configProperties).parseLoci(propPrefix);
	}

	public AlleleClassTarget[] getTargets() {
		return targets;
	}

	public static class AlleleClassLocusConfig extends LocusConfig {

		public AlleleClassLocusConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public AlleleClassLocus[] parseLoci (String propPrefix) throws AnalysisException {
			
			Locus[] baseLoci = super.parseLoci(propPrefix);
			
			AlleleClassLocus[] loci = new AlleleClassLocus[baseLoci.length];
			for (int lIdx = 0; lIdx < loci.length; lIdx++) {
				String propLocusPrefix = propPrefix+"locus."+baseLoci[lIdx].getName();
				String locusChromosome = baseLoci[lIdx].getChromosome();
				AlleleClassTarget[] targets = AlleleClassTarget.parseTargetConfig(configProperties, propLocusPrefix, locusChromosome);
				Locus baseLocus = baseLoci[lIdx];
				loci[lIdx] = new AlleleClassLocus(baseLocus.getName(), baseLocus.getReadSearchIntervalCoords(), baseLocus.getAnchors(), targets);
			}
			return loci;
		}
	}
}
