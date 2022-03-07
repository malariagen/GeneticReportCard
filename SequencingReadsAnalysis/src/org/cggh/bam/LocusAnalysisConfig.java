package org.cggh.bam;

import java.io.File;
import org.cggh.common.exceptions.AnalysisException;

public abstract class LocusAnalysisConfig extends BaseAnalysisConfig {

	public static final int DEFAULT_MAX_READ_MISMATCHES = 10;

	public static final String PROP_MAX_READ_MISMATCHES = "alignment.maxReadMismatches";
	public static final String PROP_MAX_INDEL_SIZE      = "alignment.maxIndelSize";
	
	protected Locus[] loci;

	protected boolean analyzeUnmappedReads;
	protected int     maxReadMismatches;
	protected int     maxIndelSize;
	protected boolean useBamAlignment;
	

	public LocusAnalysisConfig(File configFile, String propPrefix, boolean useBamAlignment) throws AnalysisException {
		super(configFile, propPrefix);
		
		maxReadMismatches  = this.getIntProperty(propPrefix+PROP_MAX_READ_MISMATCHES, DEFAULT_MAX_READ_MISMATCHES);	
		maxIndelSize       = this.getIntProperty(propPrefix+PROP_MAX_INDEL_SIZE,      0);	
		useBamAlignment = useBamAlignment;

		loci = parseLocusConfig ();
		analyzeUnmappedReads = false;
		for (int i = 0; i < loci.length; i++) {
			if (loci[i].getAnalyzeUnmappedReads()) {
				analyzeUnmappedReads = true;
				break;
			}
		}
	}
	
	public Locus[] getLoci() {
		return loci;
	}
	
	public int getMaxReadMismatches() {
		return maxReadMismatches;
	}
	
	public int getMaxIndelSize() {
		return maxIndelSize;
	}
	
	public boolean getUseBamAlignment () {
		return useBamAlignment;
	}
	
	public boolean getAnalyzeUnmappedReads() {
		return analyzeUnmappedReads;
	}
	
	public String getPrintableDisplay() {
	    return super.getPrintableDisplay() +
	         "\nanalyzeUnmappedReads = " + getAnalyzeUnmappedReads() +
		     "\nmaxReadMismatches = "    + getMaxReadMismatches() +
		     "\nmaxIndelSize = "         + getMaxIndelSize() +
             "\nuseBamAlignment = "      + getUseBamAlignment();    		
    }
	
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