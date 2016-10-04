package org.cggh.bam;

import java.util.Properties;

import org.cggh.common.config.BaseConfig;
import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.genome.GenomeRegion;

public class Locus {

	String       name;	
	String       readSearchIntervalCoords;
	GenomeRegion readSearchInterval;
	Anchor[]     anchors;

	public Locus(String name, String readSearchIntervalCoords, Anchor[] anchors) throws AnalysisException {
		this.name = name;
		this.readSearchIntervalCoords = readSearchIntervalCoords;
		this.readSearchInterval = GenomeRegion.parseRegion(readSearchIntervalCoords);
		this.anchors = anchors;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getReadSearchIntervalCoords() {
		return readSearchIntervalCoords;
	}

	public GenomeRegion getReadSearchInterval() {
		return readSearchInterval;
	}

	public void setReadSearchInterval(GenomeRegion readSearchInterval) {
		this.readSearchInterval = readSearchInterval;
	}

	public String getChromosome() {
		return readSearchInterval.getChromosome();
	}

	public Anchor[] getAnchors() {
		return anchors;
	}

	public void setAnchors(Anchor[] anchors) {
		this.anchors = anchors;
	}
	
	public static Locus[] parseLocusConfig(Properties configProperties, String propPrefix) throws AnalysisException {
		return new LocusConfig(configProperties).parseLoci(propPrefix);
	}

	
	public static class LocusConfig extends BaseConfig {

		public LocusConfig(Properties configProperties) {
			super(configProperties);
		}
		
		public Locus[] parseLoci (String propPrefix) throws AnalysisException {
			String[] locusNames = getStringListProperty(propPrefix+"loci");
			Locus[] loci = new Locus[locusNames.length];
			for (int i = 0; i < loci.length; i++) {
				String locusPrefix = propPrefix+"locus."+locusNames[i];
				
				String regionStr = getProperty(locusPrefix+".region");
				GenomeRegion region = GenomeRegion.parseRegion(regionStr);
				String chr = region.getChromosome();
							
				String[] anchorsStr = getStringListProperty(locusPrefix+".anchors");
				Anchor[] anchors = new Anchor[anchorsStr.length];
				for (int idx = 0; idx < anchors.length; idx++) {
					String[] parts = anchorsStr[idx].split("@");
					anchors[idx] = new Anchor(chr+":"+parts[0], parts[1]);
				}
				loci[i] = new Locus(locusNames[i], regionStr, anchors);
			}
			return loci;
		}
	}
}
