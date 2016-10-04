package org.cggh.common.sequence;

import org.cggh.common.util.TextUtilities;

public class SequenceUtilities {
	
	
	public static String makeFastaAlignment (Sequence[] sequences) {
		return makeFastaAlignment (sequences, 0, sequences.length);
	}
	
	public static String makeFastaAlignment (Sequence[] sequences, int startIdx, int count) {
	    StringBuffer sb = new StringBuffer();
	    int endIdx = count + startIdx;
        for (int i = startIdx; i < endIdx; i++) {
        	Sequence sequence = sequences[i];
        	appendFastaEntry (sb, sequence.getId(), sequence.getData());
            sb.append("\r\n");
        }       
		return sb.toString();
	}
	
	public static String makeFastaEntry (Sequence sequence) {
		return  makeFastaEntry (sequence.getId(), sequence.getData());
	}
	
	public static String makeFastaEntry (String caption, String sequence) {
	    StringBuffer sb = new StringBuffer();
	    appendFastaEntry (sb, caption, sequence);
		return sb.toString();
	}
	
	public static void appendFastaEntry (StringBuffer sb, String caption, String sequence) {

		// Replace leading and trailing spaces with gaps
		String rawSeq = sequence.replace(' ', '-'); 
        String[] seqLines = TextUtilities.split (rawSeq, 60);
        sb.append(">");
        sb.append(caption);
        for (int j = 0; j < seqLines.length; j++) {
            sb.append("\r\n");
	        sb.append(seqLines[j]);
        }
        sb.append("\r\n");        
	}
	
	public static String makeTabbedAlignmentEntry (String caption, String sequence) {
		return (caption + "\t" + sequence);
	}

	
	// *****************************************************
	
	final private static String[][] aaCodons = {
		{"AAA", "AAG"              },		// K
		{              "AAC", "AAT"},		// N
		//
		{"ACA", "ACG", "ACC", "ACT"},		// T
		{"ATA",        "ATC", "ATT"},		// I
		{       "ATG"              },		// M
		{"GAA", "GAG"              },		// E
		{              "GAC", "GAT"},		// D
		{"GGA", "GGG", "GGC", "GGT"},		// G
		{"GCA", "GCG", "GCC", "GCT"},		// A
		{"GTA", "GTG", "GTC", "GTT"},		// V
		{"CAA", "CAG"              },		// Q
		{              "CAC", "CAT"}, 		// H
		{"CGA", "CGG", "CGC", "CGT", 
		 "AGA", "AGG"              }, 		// R
		{"CCA", "CCG", "CCC", "CCT"},		// P
		{"CTA", "CTG", "CTC", "CTT", 
		 "TTA", "TTG"              },		// L
//       "TAA", "TAG"						- stop
		{              "TAC", "TAT"},		// Y
//       "TGA"								- stop
		{       "TGG"              },		// W
		{              "TGC", "TGT"},		// C
		{"TCA", "TCG", "TCC", "TCT", 
		               "AGC", "AGT"},		// S
		{              "TTC", "TTT"},		// F
	};
	
	
	/**
	 * Get the set of codons associated with a particular amino acid.
	 * 
	 * @param aa The amino acid whose codons are to be found
	 * @return an array containing the codons (or null if aa was not a valid amino acid)
	 */
	public static String[] getAaCodons (char aa) {
		switch (aa) {
		case 'k':
		case 'K':
			return aaCodons[0];
		case 'n':
		case 'N':
			return aaCodons[1];
		case 't':
		case 'T':
			return aaCodons[2];
		case 'i':
		case 'I':
			return aaCodons[3];
		case 'm':
		case 'M':
			return aaCodons[4];
		case 'e':
		case 'E':
			return aaCodons[5];
		case 'd':
		case 'D':
			return aaCodons[6];
		case 'g':
		case 'G':
			return aaCodons[7];
		case 'a':
		case 'A':
			return aaCodons[8];
		case 'v':
		case 'V':
			return aaCodons[9];
		case 'q':
		case 'Q':
			return aaCodons[10];
		case 'h':
		case 'H':
			return aaCodons[11];
		case 'r':
		case 'R':
			return aaCodons[12];
		case 'p':
		case 'P':
			return aaCodons[13];
		case 'l':
		case 'L':
			return aaCodons[14];
		case 'y':
		case 'Y':
			return aaCodons[15];
		case 'w':
		case 'W':
			return aaCodons[16];
		case 'c':
		case 'C':
			return aaCodons[17];
		case 's':
		case 'S':
			return aaCodons[18];
		case 'f':
		case 'F':
			return aaCodons[19];
		}
		return null;
	}
	
	
	public static String  getReverseComplementSequence(String sequence) {
		int len = sequence.length();
		char[] result = new char[len];
		for (int i = 0; i < len; i++) {
			char nt = sequence.charAt(len-i-1);
			result[i] = ((nt == 'A') ? 'T' : ((nt == 'T') ? 'A' : ((nt == 'C') ? 'G' : 'C')));
		}
		return new String(result);
	}
	

	public static String  translateNtSequence(String sequence) {
		StringBuffer sb = new StringBuffer(sequence.length() / 3);
		for (int idx = 0; idx < sequence.length(); idx += 3) {
			char aa = Translation.getCodon(sequence, idx+1, false).getAmino();
			sb.append(aa);
		}
		return sb.toString();
	}

}
