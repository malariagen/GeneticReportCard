package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.sampleClass.AlleleClassTarget.ClassAllele;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.AnalysisException;
import java.io.*;

public class AlleleClassAnalysisConfig extends TargetAnalysisConfig {

	public AlleleClassAnalysisConfig(File configFile, String propPrefix, boolean useAlignment) throws AnalysisException {
		super(configFile, propPrefix, useAlignment);
	}
	
	public AlleleClassLocus[] getLoci () {
		return (AlleleClassLocus[]) loci;
	}
	
	@Override
	public AlleleClassLocus[] parseLocusConfig () throws AnalysisException {
		Locus[] baseLoci = super.parseLocusConfig();
		AlleleClassLocus[] loci = new AlleleClassLocus[baseLoci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus locus = baseLoci[lIdx];
			String propLocusPrefix = propPrefix+"locus."+baseLoci[lIdx].getName();
			// If we're using the BAM alignment, then the target will be aligned on the reference chromosome.
			// If we do not use the alignment, just use a ficticious name (the name of the locus)
			String targetChromosomeName = useBamAlignment ? locus.getReadSearchInterval().getChromosome() : locus.getName();
			AlleleClassTarget[] targets = parseTargets(propLocusPrefix, targetChromosomeName);
			loci[lIdx] = new AlleleClassLocus(locus.getName(), locus.getReadSearchIntervalCoords(), locus.getAnchors(), locus.getAnalyzeUnmappedReads(), targets);
		}
		return loci;
	}
	
	@Override
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