package org.cggh.bam.grc;

import org.cggh.bam.target.*;
import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public class GrcConfig extends BaseConfig {
	
	public static final String PROP_PREFIX = "grc.";
	
	private TargetLocus[] loci;
	
	public GrcConfig (File configFile) throws AnalysisException  {
		super(configFile);
		loci = TargetLocus.parseLocusConfig(this.configProperties, PROP_PREFIX);
	}
	
	public TargetLocus[] getLoci() {
		return loci;
	}
}
