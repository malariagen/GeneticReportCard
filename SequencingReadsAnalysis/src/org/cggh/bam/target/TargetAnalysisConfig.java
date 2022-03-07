package org.cggh.bam.target;

import java.io.File;

import org.cggh.bam.*;
import org.cggh.common.exceptions.AnalysisException;

public class TargetAnalysisConfig extends LocusAnalysisConfig {

	public TargetAnalysisConfig(File configFile, String propPrefix, boolean useBamAlignment) throws AnalysisException {
		super(configFile, propPrefix, useBamAlignment);
	}
	
	public TargetLocus[] getLoci () {
		return (TargetLocus[]) loci;
	}
	
	@Override
	public TargetLocus[] parseLocusConfig () throws AnalysisException {
		Locus[] baseLoci = super.parseLocusConfig();
		TargetLocus[] loci = new TargetLocus[baseLoci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			Locus bLocus = baseLoci[lIdx];
			String propLocusPrefix = propPrefix+"locus."+bLocus.getName();			
			// If we're using the BAM alignment, then the target will be aligned on the reference chromosome.
			// If we do not use the alignment, just use a ficticious name (the name of the locus)
			String targetChromosomeName = useBamAlignment ? bLocus.getReadSearchInterval().getChromosome() : bLocus.getName();
			TargetLocus locus = new TargetLocus(bLocus.getName(), bLocus.getReadSearchIntervalCoords(), bLocus.getAnchors(), bLocus.getAnalyzeUnmappedReads());
			Target[] targets = parseTargets(propLocusPrefix, targetChromosomeName);
			locus.setTargets(targets);
			loci[lIdx] = locus;
		}
		return loci;
	}
	
	public Target[] parseTargets (String propLocusPrefix, String targetChromosomeName) throws AnalysisException {
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
				coordStrs[j] = targetChromosomeName+":"+coordStrs[j];
			}
			targets[idx] = createTarget(parts[0], coordStrs, isReverse);				
		}
		return targets;
	}
	
	protected Target createTarget (String name, String[] targetCoords, boolean isReverse) throws AnalysisException {
		return new Target (name, targetCoords, isReverse);
	}
    
	public String getPrintableDisplay() {
	    return super.getPrintableDisplay();    		
    }
}