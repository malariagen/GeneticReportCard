package org.cggh.bam.sampleClass;

import org.cggh.bam.target.*;
import org.cggh.common.exceptions.*;


public class ClassTarget extends Target {
	
	private ClassAllele[] alleles;
	
	public ClassTarget (String name, String[] targetCoords, boolean isReverse, ClassAllele[] alleles) throws AnalysisException {
		super (name, targetCoords, isReverse);
		this.alleles = alleles;
	}

	public ClassAllele[] getAlleles() {
		return alleles;
	}

	public static class ClassAllele {
		private String name;
		private String[] sequences;
		
		public ClassAllele(String name, String[] sequences) {
			this.name = name;
			this.sequences = sequences;
		}

		public String getName() {
			return name;
		}

		public String[] getSequences() {
			return sequences;
		}
	}
}


