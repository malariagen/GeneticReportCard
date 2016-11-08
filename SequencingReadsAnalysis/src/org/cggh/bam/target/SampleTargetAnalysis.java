package org.cggh.bam.target;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public abstract class SampleTargetAnalysis extends SampleLocusAnalysis {
	
	protected Target[]       allTargets;
	protected String[]       allTargetNames;
	protected HashMap<String,Integer> targetIdxTable = new HashMap<String,Integer>();

	public SampleTargetAnalysis (File refFastaFile, File chrMapFile, File outRootFolder) throws AnalysisException  {
		super (refFastaFile, chrMapFile, outRootFolder);
	}
	
	protected void registerLoci (Locus[] loci) {
		super.registerLoci (loci);
		
		// Get all the targets and their names
		ArrayList<Target> targetList = new ArrayList<Target>();
		ArrayList<String> targetNameList = new ArrayList<String>();
		
		for (int lIdx = 0; lIdx < loci.length; lIdx++) {
			TargetLocus locus = (TargetLocus)loci[lIdx];
			Target[] targets = locus.getTargets();
			for (int tIdx = 0; tIdx < targets.length; tIdx++) {
				Target target = targets[tIdx];
				String targetName = locus.getName()+"_"+target.getName();
				int allTargetIdx = targetList.size();
				targetList.add(target);
				targetNameList.add(targetName);
				targetIdxTable.put(targetName, allTargetIdx);
			}
		}
		allTargets = targetList.toArray(new Target[targetList.size()]);
		allTargetNames = targetNameList.toArray(new String[targetNameList.size()]);
	}
	
	protected int getTargetIndex (String targetName) {
		Integer idxObj = targetIdxTable.get(targetName);
		if (idxObj == null) {
			return -1;
		}
		return idxObj.intValue();
	}

}
