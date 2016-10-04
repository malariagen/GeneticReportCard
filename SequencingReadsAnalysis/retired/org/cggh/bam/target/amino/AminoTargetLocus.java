package org.cggh.bam.target.amino;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import java.util.*;


public class AminoTargetLocus extends Locus {

	private AminoTarget[] targets;
	
	public AminoTargetLocus(String name, String readSearchIntervalCoords, Anchor[] anchors, AminoTarget[] targets) throws AnalysisException {
		super (name, readSearchIntervalCoords, anchors);
		this.targets = targets;
	}

	public static AminoTargetLocus[] parseLocusConfig(Properties configProperties, String propPrefix) throws AnalysisException {
		return new NtAlleleTargetLocusConfig(configProperties).parseLoci(propPrefix);
	}

	public AminoTarget[] getTargets() {
		return targets;
	}

	public static class NtAlleleTargetLocusConfig extends LocusConfig {

		public NtAlleleTargetLocusConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public AminoTargetLocus[] parseLoci (String propPrefix) throws AnalysisException {
			
			Locus[] baseLoci = super.parseLoci(propPrefix);
			
			AminoTargetLocus[] loci = new AminoTargetLocus[baseLoci.length];
			for (int lIdx = 0; lIdx < loci.length; lIdx++) {
				String propLocusPrefix = propPrefix+"locus."+baseLoci[lIdx].getName();
				String locusChromosome = baseLoci[lIdx].getChromosome();
				AminoTarget[] targets = AminoTarget.parseTargetConfig(configProperties, propLocusPrefix, locusChromosome);
				Locus baseLocus = baseLoci[lIdx];
				loci[lIdx] = new AminoTargetLocus(baseLocus.getName(), baseLocus.getReadSearchIntervalCoords(), baseLocus.getAnchors(), targets);
			}
			return loci;
		}
	}
	
	
	
}
