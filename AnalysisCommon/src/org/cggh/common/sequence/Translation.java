package org.cggh.common.sequence;

import java.io.*;

import org.cggh.common.exceptions.AnalysisException;


public class Translation {
	
	// Translation table
	// NT1: AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
	// NT2: AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
	// NT3: ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
	// AA:  KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV.Y.YSSSS.CWCLFLF      
	
	public static final Codon[] codons = loadCodons();
	
	public static char getComplementaryNt (char nt) {
		switch (nt) {
		case 'A':	return 'T';
		case 'C':	return 'G';
		case 'G':	return 'C';
		case 'T':	return 'A';
		}
		return 'X';
	}

	public static Codon getCodon (String codon) {
		int idx = getCodonIndex (codon);
		if (idx < 0) {
			return null;
		}
		return codons[idx];
	}
	
	public static Codon getCodon (char nt1, char nt2, char nt3) {
		int idx = getCodonIndex (nt1, nt2, nt3);
		if (idx < 0) {
			return null;
		}
		return codons[idx];
	}
	
	public static Codon getCodon (String sequence, int startPos, boolean isReverseStrand) {
		int startPosIdx = startPos-1;
		int idx = isReverseStrand ?
		        getCodonIndex (getComplementaryNt(sequence.charAt(startPosIdx-0)), 
		        		       getComplementaryNt(sequence.charAt(startPosIdx-1)), 
		        			   getComplementaryNt(sequence.charAt(startPosIdx-2))) :
				getCodonIndex (sequence.charAt(startPosIdx+0), 
						       sequence.charAt(startPosIdx+1), 
						       sequence.charAt(startPosIdx+2));
		if (idx < 0) {
			return null;
		}
		return codons[idx];
	}
	
	static int getCodonIndex (String codon) {
		return getCodonIndex (codon.charAt(0), codon.charAt(1), codon.charAt(2));
	}
	
	static int getCodonIndex (char nt1, char nt2, char nt3) {
		int idx1 = getNtIndex(nt1);
		int idx2 = getNtIndex(nt2);
		int idx3 = getNtIndex(nt3);
		if ((idx1 < 0) || (idx2 < 0) || (idx3 < 0)) {
			return -1;
		}
		int idx = (16 * idx1) + (4 * idx2) + idx3;
		return idx;
	}
	
	
	static int getNtIndex (char nt) {
		char n = Character.toUpperCase(nt);
		return ((n == 'A') ? 0 : (n == 'C' ? 1 : (n == 'G' ? 2 : (n == 'T' ? 3 : -1))));
	}
	
	static final String NT = "ACGT";
	static final String AA_TRANS  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV.Y.YSSSS.CWCLFLF";
	static final String SYN_4_WAY = "0000111100000000000011111111111100001111111111110000111100000000";

	private static Codon[] loadCodons() {
		
		// Create the basic list of codons
		Codon[] allCodons = new Codon[64]; 
		for (int p1 = 0; p1 < 4; p1++) {
			for (int p2 = 0; p2 < 4; p2++) {
				for (int p3 = 0; p3 < 4; p3++) {
					int idx = (16 * p1) + (4 * p2) + p3;
					char aa = AA_TRANS.charAt(idx);
					boolean fourWaySynonymous = (SYN_4_WAY.charAt(idx) == '1');
					String cod = "" + NT.charAt(p1) + NT.charAt(p2) + NT.charAt(p3);
					Codon codon = new Codon (cod, aa, fourWaySynonymous);
					allCodons[idx] = codon;
				}
			}
		}
		//for (int codonIdx = 0; codonIdx < 64; codonIdx++) {
		//	System.out.println(allCodons[codonIdx]);
		//}
		
		// For each codon, work out all the possible mutations
		for (int codonIdx = 0; codonIdx < 64; codonIdx++) {
			Codon codon = allCodons[codonIdx];
			Mutation[][] codMutations = codon.getMutationsTable();
			for (int pos = 0; pos < 3; pos++) {
				int posShiftFactor = (pos == 0) ? 16 : ((pos == 1) ? 4 : 1);
				int posMask = 63 - (3 * posShiftFactor);
				for (int ntIdx = 0; ntIdx < 4; ntIdx++) {
					int a = (codonIdx & posMask);
					int b = (ntIdx * posShiftFactor);
					int mutCodonIdx = a + b;
					//int mutCodonIdx = (codonIdx & posMask) + (ntIdx * posFactor);
					Codon mutCodon = allCodons[mutCodonIdx];
					codMutations[pos][ntIdx] = new Mutation(codon, mutCodon, pos);
				}
			}
		}
		return allCodons;
	}
	
	
	/**
	 * 
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		try {
			File outfile = (args.length == 0) ? null : new File (args[0]);
			PrintWriter w = (outfile == null) ? new PrintWriter(System.out) : new PrintWriter (outfile);
			
			StringBuffer sb = new StringBuffer(10000);
			for (int codonIdx = 0; codonIdx < 64; codonIdx++) {
				Codon codon = codons[codonIdx];
				sb.append("#");
				sb.append(codonIdx);
				sb.append(": ");
				sb.append(codon);
				for (int pos = 0; pos < 3; pos++) {
					for (int ntIdx = 0; ntIdx < 4; ntIdx++) {
						// Test the getMutation() method
						char nt = NT.charAt(ntIdx);
						Mutation mut = codon.getMutation(nt, pos);
						if (mut.getWildType() != codon) {
							throw new AnalysisException ("Bad mutation found: " + mut + " for " + codon);
						}
						sb.append("\n  ");
						sb.append(mut.toString());
					}
				}
				w.println(sb.toString());
				sb.setLength(0);
			}
			w.close();
			
		} catch (Exception e) {
			System.err.println("Error writing codon table out: "+e);
		}
	}

}