package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.AnalysisException;
import java.io.*;

public class AlleleClassAnalysisConfig extends TargetAnalysisConfig {

	public AlleleClassAnalysisConfig(File configFile, String propPrefix) throws AnalysisException {
		super(configFile, propPrefix);
	}
	
	public AlleleClassLocus[] getLoci () {
		return (AlleleClassLocus[]) loci;
	}
	
	@Override
	public AlleleClassLocus[] parseLocusConfig () throws AnalysisException {
		Locus[] baseLoci = super.parseLocusConfig();
		AlleleClassLocus[] loci = new AlleleClassLocus[baseLoci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			String propLocusPrefix = propPrefix+"locus."+baseLoci[lIdx].getName();
			String locusChromosome = baseLoci[lIdx].getChromosome();
			AlleleClassTarget[] targets = AlleleClassTarget.parseTargetConfig(configProperties, propLocusPrefix, locusChromosome);
			Locus baseLocus = baseLoci[lIdx];
			loci[lIdx] = new AlleleClassLocus(baseLocus.getName(), baseLocus.getReadSearchIntervalCoords(), baseLocus.getAnchors(), baseLocus.getAnalyzeUnmappedReads(), targets);
		}
		return loci;
	}
}