package org.cggh.common.counters;

public class AminoAlleleCounter extends AlleleCounter {

	private static final char[] ALLELES = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','.'};
	
	public AminoAlleleCounter () { 
		super(ALLELES);
	}
	
	protected int getIndex(char allele) {
		switch(allele) {
		case 'A':	return 0;
		case 'C':	return 1;
		case 'D':	return 2;
		case 'E':	return 3;
		case 'F':	return 4;
		case 'G':	return 5;
		case 'H':	return 6;
		case 'I':	return 7;
		case 'K':	return 8;
		case 'L':	return 9;
		case 'M':	return 10;
		case 'N':	return 11;
		case 'P':	return 12;
		case 'Q':	return 13;
		case 'R':	return 14;
		case 'S':	return 15;
		case 'T':	return 16;
		case 'V':	return 17;
		case 'W':	return 18;
		case 'Y':	return 19;
		case '.':	return 20;
		}
		return -1;
	}
}
