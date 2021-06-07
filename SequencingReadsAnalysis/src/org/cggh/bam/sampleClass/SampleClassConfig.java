package org.cggh.bam.sampleClass;

import org.cggh.bam.*;
import org.cggh.bam.sampleClass.ClassTarget.ClassAllele;
import org.cggh.bam.target.*;
import org.cggh.common.exceptions.AnalysisException;
import java.io.*;

public class SampleClassConfig extends TargetAnalysisConfig {
	
	public static final String PROP_PREFIX = "sampleClass.";
	
	protected String[]           classes;

	public SampleClassConfig(File configFile) throws AnalysisException {
		super(configFile, PROP_PREFIX, false);
		classes = getStringListProperty(propPrefix+"classes");
	}
	
	public ClassLocus[] getLoci () {
		return (ClassLocus[]) loci;
	}
	
	public String[] getClasses() {
		return classes;
	}

	@Override
	public ClassLocus[] parseLocusConfig () throws AnalysisException {
		Locus[] baseLoci = super.parseLocusConfig();
		ClassLocus[] loci = new ClassLocus[baseLoci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus bLocus = baseLoci[lIdx];
			String propLocusPrefix = propPrefix+"locus."+bLocus.getName();
			// If we're using the BAM alignment, then the target will be aligned on the reference chromosome.
			// If we do not use the alignment, just use a ficticious name (the name of the locus)
			String targetChromosomeName = useBamAlignment ? bLocus.getReadSearchInterval().getChromosome() : bLocus.getName();
			ClassLocus locus = new ClassLocus(bLocus.getName(), bLocus.getReadSearchIntervalCoords(), bLocus.getAnchors(), bLocus.getAnalyzeUnmappedReads());
			ClassTarget[] targets = parseTargets(propLocusPrefix, targetChromosomeName);
			locus.setTargets(targets);
			loci[lIdx] = locus;
		}
		return loci;
	}
	
	@Override
	public ClassTarget[] parseTargets (String propLocusPrefix, String locusChromosome) throws AnalysisException {
		Target[] baseTargets = super.parseTargets(propLocusPrefix, locusChromosome);
		
		ClassTarget[] targets = new ClassTarget[baseTargets.length];
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
			targets[tIdx] = new ClassTarget(baseTarget.getName(), baseTarget.getTargetCoords(), baseTarget.isReverse(), alleles);
		}
		return targets;
		
	}
}