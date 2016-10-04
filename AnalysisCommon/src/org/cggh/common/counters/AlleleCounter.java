
package org.cggh.common.counters;

import java.util.Arrays;

public abstract class AlleleCounter {
	
	protected AlleleCount[] counts;
	
	public AlleleCounter (char[] alleles) {
		counts = new AlleleCount[alleles.length];
		for (int i = 0; i < counts.length; i++) {
			counts[i] = new AlleleCount(alleles[i]);
		}
	}
	
	protected abstract int getIndex (char allele);
	
	public void reset () {
		for (int i = 0; i < counts.length; i++) {
			counts[i].count = 0;
		}
	}
	
	public boolean isValidAllele(char allele) {
		return (getIndex (allele) >= 0);
	}

	public int getCount(char allele) {
		int idx = getIndex(allele);
		return counts[idx].count;
	}
	
	public void setCount(char allele, int newCount) {
		int idx = getIndex(allele);
		counts[idx].count = newCount;
	}
	
	public void setCounts(int[] newCounts) {
		for (int i = 0; i < counts.length; i++) {
			counts[i].count = newCounts[i];
		}
	}
	
	public void addCount(char allele, int incrCount) {
		int idx = getIndex(allele);
		counts[idx].count += incrCount;
	}
	
	public void increment(char allele) {
		int idx = getIndex(allele);
		counts[idx].count++;
	}
	
	/**
	 * Returns the total number of reads foor this sample at this SNP (i.e. the total coverage.
	 * 
	 * @return
	 */
	public int getCumulativeCount() {
		int cumulativeCount = 0;
		for (int i = 0; i < counts.length; i++) {
			cumulativeCount += counts[i].count;
		}
		return cumulativeCount;
	}

	/**
	 * Returns the alleles, in decreasing order by number of reads.
	 * The order is determined by getSortedAlleleCounts()
	 * 
	 * @return the array of AlleleCount objects, one for each allele
	 */
	public char[] getSortedAlleles() {
		AlleleCount[] sortedCounts = getSortedAlleleCounts();
		char[] result = new char[sortedCounts.length];
		for (int i = 0; i < sortedCounts.length; i++) {
			result[i] = sortedCounts[i].allele;
		}
		return result;
	}

	/**
	 * Returns the alleles counters, in decreasing order by number of reads.
	 * The order is determined by getSortedAlleleCounts()
	 * The order for alleles with the same number of reads is undefined,
	 * but specific behaviours may be implemented by subclasses.
	 * 
	 * @return the array of AlleleCount objects, one for each allele
	 */
	public AlleleCount[] getSortedAlleleCounts() {
		AlleleCount[] sortedCounts = getAlleleCounts();
		Arrays.sort(sortedCounts);
		return sortedCounts;
	}
	
	public AlleleCount[] getAlleleCounts() {
		AlleleCount[] alleleCounts = new AlleleCount[counts.length];
		for (int i = 0; i < counts.length; i++) {
			alleleCounts[i] = counts[i];
		}
		return alleleCounts;
	}
	
	
	public class AlleleCount implements Comparable<AlleleCount> {
		
		char allele;
		int  count;
		
		public AlleleCount (char allele) {
			this (allele, 0);
		}
		
		public AlleleCount (char allele, int count) {
			this.allele = allele;
			this.count = count;
		}
		
		public char getAllele () {
			return allele;
		}

		public int getCount () {
			return count;
		}

		@Override
		public int compareTo(AlleleCount o) {
			return o.count - count;  // Reverse order
		}
	}
}
