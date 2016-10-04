package org.cggh.common.sequence;


public class Codon {
	
	private String       sequence;
	private char         amino;
	private Mutation[][] mutations = new Mutation[3][4];
	private boolean      fourWaySynonymous;
	
	public Codon(String sequence, char amino, boolean fourWaySynonymous) {
		this.sequence = sequence;
		this.amino = amino;
		this.fourWaySynonymous = fourWaySynonymous;
	}
	
	public String getSequence() {
		return sequence;
	}

	public char getNt(int pos) {
		int posIdx = pos-1;
		return sequence.charAt(posIdx);
	}

	public char getAmino() {
		return amino;
	}

	public Mutation getMutation(char nt, int pos) {
		int posIdx = pos-1;
		int ntIdx = Translation.getNtIndex(nt);
		return mutations[posIdx][ntIdx];
	}
	
	public boolean isStopCodon() {
		return (amino == '.');
	}
	
	public boolean isFourWaySynonymous() {
		return fourWaySynonymous;
	}
	
	public String toString() {
		return ("Codon: " + sequence + " (" + (isStopCodon() ? "STOP" : amino) + ")");
	}

	Mutation[][] getMutationsTable() {
		return mutations;
	}
}

