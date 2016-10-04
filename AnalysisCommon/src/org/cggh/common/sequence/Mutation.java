package org.cggh.common.sequence;


public class Mutation {
	
	private Codon wildType;
	private Codon mutated;
	private int   pos;     // 0-based indexing
	
	public Mutation(Codon wildType, Codon mutated, int pos) {
		this.wildType = wildType;
		this.mutated = mutated;
		this.pos = pos;
	}
	
	public Codon getWildType() {
		return wildType;
	}

	public Codon getMutated() {
		return mutated;
	}

	public int getPos() {
		return pos;
	}

	public char getWildAllele() {
		return wildType.getNt(pos);
	}

	public char getMutatedAllele() {
		return mutated.getNt(pos);
	}

	public String getWildSequence() {
		return wildType.getSequence();
	}

	public String getMutatedSequence() {
		return mutated.getSequence();
	}

	public char getWildAmino() {
		return wildType.getAmino();
	}

	public char getMutatedAmino() {
		return mutated.getAmino();
	}

	public boolean isSynonymous() {
		return (getWildAmino() == getMutatedAmino());
	}
	
	public String toString() {
		return ("Mutation: pos=" + (pos+1) + " " 
				+ getWildAllele() + "->" + getMutatedAllele() 
				+ " (" + getWildSequence() + "->" +  getMutatedSequence() + ", "
				+ (wildType.isStopCodon()? "STOP" : getWildAmino()) + "->" + (mutated.isStopCodon()? "STOP" : getMutatedAmino()) + ") ");
	}
}
