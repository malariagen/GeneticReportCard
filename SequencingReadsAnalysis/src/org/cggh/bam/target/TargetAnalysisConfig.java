package org.cggh.bam.target;

import java.io.File;

import org.cggh.bam.Locus;
import org.cggh.bam.LocusAnalysisConfig;
import org.cggh.common.exceptions.AnalysisException;

public class TargetAnalysisConfig extends LocusAnalysisConfig {

	public TargetAnalysisConfig(File configFile, String propPrefix) throws AnalysisException {
		super(configFile, propPrefix);
	}
	
	public TargetLocus[] getLoci () {
		return (TargetLocus[]) loci;
	}
	
	@Override
	public TargetLocus[] parseLocusConfig () throws AnalysisException {
		Locus[] baseLoci = super.parseLocusConfig();
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