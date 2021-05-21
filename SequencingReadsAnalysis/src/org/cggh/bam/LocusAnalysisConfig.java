package org.cggh.bam;

import java.io.File;
import org.cggh.common.exceptions.AnalysisException;

public class LocusAnalysisConfig extends BaseAnalysisConfig {

	public LocusAnalysisConfig(File configFile, String propPrefix, boolean useBamAlignment) throws AnalysisException {
		super(configFile, propPrefix, useBamAlignment);
	}
	
	@Override
	public Locus[] parseLocusConfig () throws AnalysisException {
		String[] locusNames = getStringListProperty(propPrefix+"loci");
		Locus[] loci = new Locus[locusNames.length];
		for (int i = 0; i < loci.length; i++) {
			String locusPrefix = propPrefix+"locus."+locusNames[i];
			
			String regionStr = getProperty(locusPrefix+".region");
			//GenomeRegion region = GenomeRegion.parseRegion(regionStr);
			//String chr = region.getChromosome();

			String[] anchorsStr = getStringListProperty(locusPrefix+".anchors");
			Anchor[] anchors = new Anchor[anchorsStr.length];
			for (int idx = 0; idx < anchors.length; idx++) {
				String[] parts = anchorsStr[idx].split("@");
				anchors[idx] = new Anchor(Integer.parseInt(parts[0]), parts[1]);
			}
			boolean analyzeUnmappedReads = this.getBooleanProperty(locusPrefix+".analyzeUnmappedReads", false);
			loci[i] = new Locus(locusNames[i], regionStr, anchors, analyzeUnmappedReads);
		}
		return loci;
	}
}