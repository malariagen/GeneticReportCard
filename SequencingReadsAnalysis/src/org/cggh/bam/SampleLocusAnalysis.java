package org.cggh.bam;

import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.*;

public abstract class SampleLocusAnalysis extends SampleAnalysis {

	protected Locus[]  loci;
	protected String[] locusNames;
	protected HashMap<String,Integer> locusIdxTable = new HashMap<String,Integer>();
	
	public SampleLocusAnalysis (File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, chrMapFile, outRootFolder);		
	}

	protected void registerLoci (Locus[] loci) {
		this.loci = loci;
		this.locusNames = new String[loci.length];
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			String locusName = loci[lIdx].getName();
			locusNames[lIdx] = locusName;
			locusIdxTable.put(locusName, lIdx);
		}
	}
	
	protected int getLocusIndex (String locusName) {
		Integer idxObj = locusIdxTable.get(locusName);
		if (idxObj == null) {
			return -1;
		}
		return idxObj.intValue();
	}
}
