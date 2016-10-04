package org.cggh.bam.species;

import org.cggh.bam.target.alleleClasses.AlleleClassLocus;
import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.*;

public class SpeciesConfig extends BaseConfig {
	
	public static final String PROP_PREFIX = "species.";
	
	private AlleleClassLocus[] loci;
	private String[] classes;
	
	
	public SpeciesConfig (File configFile) throws AnalysisException  {
		super(configFile);
		loci = AlleleClassLocus.parseLocusConfig(this.configProperties, PROP_PREFIX);
		classes = getStringListProperty(PROP_PREFIX+"classes");
	}
	
	public AlleleClassLocus[] getLoci() {
		return loci;
	}
	
	public String[] getClasses() {
		return classes;
	}
}
