package org.cggh.common.counters;


public class ExtendedNtAlleleCounter extends AlleleCounter {

	public static final char[] ALLELES = {'A','C','G','T','*'};
	
	public ExtendedNtAlleleCounter () { 
		super(ALLELES);
	}
	
	public ExtendedNtAlleleCounter (AlleleCount[] initCounters) {
		super(ALLELES);
		initialize (initCounters);
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
	
	public void setCounts(int a, int c, int g, int t, int spanDel) {
		counts[0].count = a;
		counts[1].count = c;
		counts[2].count = g;
		counts[3].count = t;
		counts[4].count = spanDel;
	}

	public void addCounts(int a, int c, int g, int t, int spanDel) {
		counts[0].count += a;
		counts[1].count += c;
		counts[2].count += g;
		counts[3].count += t;
		counts[4].count += spanDel;
	}
	
	public int[] getCounts() {
		return new int[] {counts[0].count, counts[1].count, counts[2].count, counts[3].count, counts[4].count};
	}
	

}
