package org.cggh.common.counters;


public class NtAlleleCounter extends AlleleCounter {

	private static final char[] ALLELES = {'A','C','G','T'};
	
	public NtAlleleCounter () { 
		super(ALLELES);
	}
	
	public NtAlleleCounter (AlleleCount[] initCounters) {
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
		}
		return -1;
	}
		
	public void setCounts(int a, int c, int g, int t) {
		counts[0].count = a;
		counts[1].count = c;
		counts[2].count = g;
		counts[3].count = t;
	}

	public void addCounts(int a, int c, int g, int t) {
		counts[0].count += a;
		counts[1].count += c;
		counts[2].count += g;
		counts[3].count += t;
	}
	
	public int[] getCounts() {
		return new int[] {counts[0].count, counts[1].count, counts[2].count, counts[3].count};
	}
	
}
