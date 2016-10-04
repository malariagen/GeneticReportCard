package org.cggh.bam.genotyping;

import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import java.io.*;

public class GenotypingConfig extends BaseConfig {
	
	private Locus[] loci;
	private int maxThreads;
	
	
	/* EXAMPLE
	#genotype.loci=pfcrt_core,pfcrt_exon4
	#
	#genotype.locus.pfcrt_core.region=Pf3D7_07_v3:403500-403800
	#genotype.locus.pfcrt_core.targets=pfcrt_72-76@403612-403626
	#genotype.locus.pfcrt_core.anchors=403593@TATTATTTATTTAAGTGTA,403627@ATTTTTGCTAAAAGAAC
	#
	#genotype.locus.pfcrt_exon4.region=Pf3D7_07_v3:404250-404550
	#genotype.locus.pfcrt_exon4.targets=pfcrt_220@404407-404409
	#genotype.locus.pfcrt_exon4.anchors=404379@TATCATATTTAATCTTGTCTTAATTAGT
	*/

	public GenotypingConfig (File configFile) throws AnalysisException  {
		super(configFile);
		
		String[] locusNames = getStringListProperty("genotype.loci");
		loci = new Locus[locusNames.length];
		for (int i = 0; i < loci.length; i++) {
			String prefix = "genotype.locus."+locusNames[i];
			
			String regionStr = getProperty(prefix+".region");
			GenomeRegion region = GenomeRegion.parseRegion(regionStr);
			String chr = region.getChromosome();
			
			String[] targetsStr = getStringListProperty(prefix+".targets");
			Target[] targets = new Target[targetsStr.length];
			for (int idx = 0; idx < targets.length; idx++) {
				String[] parts = targetsStr[idx].split("@");
				String coordStr = parts[1].trim();
				boolean isReverse = false;
				if (coordStr.startsWith("-")) {
					isReverse = true;
					coordStr = coordStr.substring(1);  // Remove -
				}
				targets[idx] = new Target(parts[0], chr+":"+coordStr, isReverse);
			}
			
			String[] anchorsStr = getStringListProperty(prefix+".anchors");
			Anchor[] anchors = new Anchor[anchorsStr.length];
			for (int idx = 0; idx < anchors.length; idx++) {
				String[] parts = anchorsStr[idx].split("@");
				anchors[idx] = new Anchor(chr+":"+parts[0], parts[1]);
			}
			loci[i] = new Locus(locusNames[i], regionStr, targets, anchors);
		}
		
		maxThreads = getIntProperty("genotype.maxThreads", 0);
	}
	
	public Locus[] getLoci() {
		return loci;
	}

	public int getMaxThreads() {
		return maxThreads;
	}

	

	/*
	 * FOR DEBUG
	 * 
	public GenotypingConfig () throws AnalysisException  {
		super (new Properties());
		locusConfigs = new Locus[] {
	        new Locus ("pfcrt_core", "Pf3D7_07_v3:403500-403800",
	    	    new Target[] {
	    	        new Target ("pfcrt_72-76", "Pf3D7_07_v3:403612-403626"),
	    	    },
			    new Anchor[] {
				    new Anchor ("Pf3D7_07_v3:403593", "TATTATTTATTTAAGTGTA"),
				    new Anchor ("Pf3D7_07_v3:403627", "ATTTTTGCTAAAAGAAC"),
			    }
		    ),
		    new Locus ("pfcrt_exon4", "Pf3D7_07_v3:404250-404550",
		        new Target[] {
			        new Target ("pfcrt_220", "Pf3D7_07_v3:404407-404409"),
			    },
			    new Anchor[] {
			        new Anchor ("Pf3D7_07_v3:404379", "TATCATATTTAATCTTGTCTTAATTAGT"),
			    }
		    ),
	    };
    }
	 */

}
