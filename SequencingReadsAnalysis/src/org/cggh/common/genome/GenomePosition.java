package org.cggh.common.genome;

import org.cggh.common.exceptions.AnalysisException;


public class GenomePosition extends GenomeRegion {

	public GenomePosition(String chr, int pos) {
		setChromosome(chr);
		setPos(pos);
	}

	public int getPos() {
		return getStartPos();
	}
	
	public void setPos(int pos) {
		setStartPos(pos);
		setStopPos(pos);
	}


	public static GenomePosition parsePosition (String gpStr) throws AnalysisException {
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
		GenomePosition p = new GenomePosition(chr, pos);
		return p;
	}
}
