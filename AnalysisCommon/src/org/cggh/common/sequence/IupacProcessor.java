package org.cggh.common.sequence;

import java.util.Arrays;

public class IupacProcessor {

	public static final char     UNKNOWN = '?';
	public static final char     NO_CALL = '-';
	
	private static final char[][] IUPAC_ALLELES;
	private static final char[]   IUPAC_CODES;
	
	
	public static char[] getEquivalentNucleotides(char iupac) {
		return IUPAC_ALLELES[iupac];
	}

	public static char getIupacCode(char[] alleles) {
		int allelesIndex = getAllelesIndex (alleles);
		if (allelesIndex < 0) {
			return UNKNOWN;
		}
		return IUPAC_CODES[allelesIndex];
	}

	public static boolean isDegenerateCode (char iupac) {
		char[] equiv = IUPAC_ALLELES[iupac];
		return ((equiv != null) && (equiv.length > 1));
	}

	static {
		
		// Create lookup tables and initialize everything to "unknown"
		char[] unknownCode = new char[] {UNKNOWN};
		IUPAC_ALLELES = new char[128][];
		Arrays.fill(IUPAC_ALLELES, unknownCode);
		
		IUPAC_CODES = new char[16];
		Arrays.fill(IUPAC_CODES, UNKNOWN);
		
		// The no-call character just translates to itself
		indexIupac (NO_CALL, new char[] {NO_CALL});
		
		// Map the IUPAC codes as appropriate
		indexIupac ('A', new char[] {'A'});
		indexIupac ('C', new char[] {'C'});
		indexIupac ('G', new char[] {'G'});
		indexIupac ('T', new char[] {'T'});
		indexIupac ('W', new char[] {'A', 'T'});
		indexIupac ('S', new char[] {'C', 'G'});
		indexIupac ('R', new char[] {'A', 'G'});
		indexIupac ('Y', new char[] {'C', 'T'});
		indexIupac ('M', new char[] {'A', 'C'});
		indexIupac ('K', new char[] {'G', 'T'});
		indexIupac ('B', new char[] {'C', 'G', 'T'});
		indexIupac ('D', new char[] {'A', 'G', 'T'});
		indexIupac ('H', new char[] {'A', 'C', 'T'});
		indexIupac ('V', new char[] {'A', 'C', 'G'});
		indexIupac ('N', new char[] {'A', 'C', 'G', 'T'});
	}
	
	private static void indexIupac (char code, char[] alleles) {
		// Sort the alleles
		Arrays.sort(alleles);
		
		// Add alleles for this code to the IUPAC_ALLELES table 
		IUPAC_ALLELES[code] = alleles;
		IUPAC_ALLELES[Character.toLowerCase(code)] = alleles;
		
		// Index the IUPAC code in the IUPAC_CODES table 
		int allelesIndex = getAllelesIndex (alleles);
		if (allelesIndex >= 0) {
			IUPAC_CODES[allelesIndex] = code;
		}
	}
	
	private static int getAllelesIndex (char[] alleles) {
		// If there are no alleles, point to NO_CALL (dash)
		if ((alleles == null) || (alleles.length == 0)) {
			return 0;
		}
		// If there is only a dash, point to NO_CALL (dash)
		if ((alleles.length == 1) && (alleles[0] == NO_CALL))  {
			return 0;
		}
		
		Arrays.sort(alleles);
		int idx = 0;
		for (int i = 0; i < alleles.length; i++) {
			char a = alleles[i];
			int alleleIdx = (a == 'A') ? 0 : (a == 'C' ? 1 : (a == 'G' ? 2 : (a == 'T' ? 3 : -1)));
			if (alleleIdx < 0) {
				return -1;  //Error
			}
			int alleleFlag = (1 << alleleIdx);
			idx |= alleleFlag;
		}
		return idx;
	}
}
