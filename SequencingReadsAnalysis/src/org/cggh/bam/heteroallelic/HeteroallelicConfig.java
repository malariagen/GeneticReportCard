
package org.cggh.bam.heteroallelic;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import java.io.*;

public class HeteroallelicConfig extends BaseAnalysisConfig {
	
    public static final String PROP_PREFIX = "heteroallelic.";
    
    private HeteroallelicLocus[] loci;

    public HeteroallelicConfig(File configFile) throws AnalysisException {
        super(configFile, PROP_PREFIX);
        
        String[] locusNames = getStringListProperty(PROP_PREFIX + "loci");
        loci = new HeteroallelicLocus[locusNames.length];
        for (int i = 0; i < loci.length; i++) {
            String locusPrefix = PROP_PREFIX + "locus." + locusNames[i];
            String regionStr = getProperty(locusPrefix + ".region");
            GenomeRegion region = GenomeRegion.parseRegion(regionStr);
            boolean isReverse = getBooleanProperty(locusPrefix + ".reverse");
            int startCodon = getIntProperty(locusPrefix + ".startCodon");
            loci[i] = new HeteroallelicLocus(locusNames[i], region, startCodon, isReverse);
        }
    }

    public HeteroallelicLocus[] getLoci() {
        return loci;
    }
    
	public String getPrintableDisplay() {
	    return super.getPrintableDisplay();    		
    }
}

