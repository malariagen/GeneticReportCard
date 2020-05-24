package org.cggh.common.genome;

import org.cggh.common.exceptions.*;


public class GenomePosition extends GenomeRegion {

	public GenomePosition(String chr, int pos) {
		setChromosome(chr);
		setPos(pos);
	}

	public GenomePosition(String gpStr) throws AnalysisException {
		GenomeCoordinates coords = parseCoordinates (gpStr);
		setChromosome(coords.chr);
		setPos(coords.pos);
	}
	
	public int getPos() {
		return getStartPos();
	}
	
	public void setPos(int pos) {
		setStartPos(pos);
		setStopPos(pos);
	}
	
    public boolean isSamePositionAs (GenomePosition other) {
		return (compareTo(other) == 0);
	}

	public static GenomePosition parsePosition (String gpStr) throws AnalysisException {
		return new GenomePosition(gpStr);
	}
	
	protected static GenomeCoordinates parseCoordinates (String gpStr) throws AnalysisException {
		// String must be of the form "<chromosome>:<pos>"
		String[] parts = gpStr.split(":");
		if (parts.length != 2) {
			throw new AnalysisException("Error parsing genomic position '" + gpStr + "': one ':' separator is required.");			
		}
		String chr = parts[0];
		String posStr = parts[1];
		int pos;
		try {
			pos = Integer.parseInt(posStr);
		} catch (Exception e) {
			throw new AnalysisException("Error parsing position of '" + gpStr + "'");			
		}
		if (pos <= 0) {
			throw new AnalysisException("Error parsing genomic position '" + gpStr + "': invalid position.");			
		}
		return new GenomeCoordinates(chr, pos);
	}
	
	static protected class GenomeCoordinates {
		String chr;
		int    pos;
		public GenomeCoordinates(String chr, int pos) {
			this.chr = chr;
			this.pos = pos;
		}
	}
}
