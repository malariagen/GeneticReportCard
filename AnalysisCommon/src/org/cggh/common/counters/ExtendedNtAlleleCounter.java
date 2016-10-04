package org.cggh.common.counters;

public class ExtendedNtAlleleCounter extends AlleleCounter {

	public static final char[] ALLELES = {'A','C','G','T','*'};
	
	public ExtendedNtAlleleCounter () { 
		super(ALLELES);
	}
	
	protected int getIndex(char allele) {
		switch(allele) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case '*':
			return 4;
		}
		return -1;
	}
}
