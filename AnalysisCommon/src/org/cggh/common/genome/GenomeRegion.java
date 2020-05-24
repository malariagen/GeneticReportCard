package org.cggh.common.genome;

import java.util.Comparator;

import org.cggh.common.exceptions.AnalysisException;

public class GenomeRegion implements Comparable<GenomeRegion> {
	
	protected String chromosome;
	protected int    startPos;
	protected int    stopPos;
	protected String name;
	
	
	public GenomeRegion () {
	}

	public GenomeRegion (String chromosome, int startPos, int stopPos) {
		setChromosome(chromosome);
		setStartPos(startPos);
		setStopPos(stopPos);
	}

	public GenomeRegion (GenomeRegion other) {
		this(other.chromosome, other.startPos, other.stopPos);
	}

	public String getChromosome() {
		return chromosome;
	}
	
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	
	public int getStartPos() {
		return startPos;
	}
	
	public void setStartPos(int startPos) {
		this.startPos = startPos;
	}
	
	public int getStopPos() {
		return stopPos;
	}
	
	public void setStopPos(int stopPos) {
		this.stopPos = stopPos;
	}

	public String getName() {
		return toString();
	}
	
	public int getSize() {
		return (1 + stopPos - startPos);
	}

	public boolean containsPos (int pos) {
		return ((pos >= startPos) && (pos <= stopPos));
	}
	
	public boolean containsPos (GenomePosition pos) {
		return (this.overlaps(pos));
	}
	
	public boolean overlaps (GenomeRegion other) {
		if (!other.chromosome.equals(this.chromosome)) {
			return false;
		}
		if (other.startPos > this.stopPos) {
			return false;
		}
		if (other.stopPos < this.startPos) {
			return false;
		}
		return true;
	}
	
	public String toString() {
		if (name == null) {
			StringBuffer sb = new StringBuffer();
			sb.append(chromosome);
			sb.append(':');
			sb.append(startPos);
			if (stopPos != startPos) {
				sb.append('-');
				sb.append(stopPos);
			}
			name = sb.toString();			
		}
		return name;
	}
	
	public int compareTo (GenomeRegion other) {
		return forwardStrainSorter.compare(this, other);
	}
	
	public static Comparator<GenomeRegion> forwardStrainSorter = new Comparator<GenomeRegion>() {
		public int compare(GenomeRegion o1, GenomeRegion o2) {
			int chrDiff = o1.chromosome.compareTo(o2.chromosome);
			if (chrDiff != 0) {
				return chrDiff;
			}
			return (o1.startPos - o2.startPos);
		}
	};

	public static Comparator<GenomeRegion> reverseStrainSorter = new Comparator<GenomeRegion>() {
		public int compare(GenomeRegion o1, GenomeRegion o2) {
			int chrDiff = o2.compareTo(o1);
			if (chrDiff != 0) {
				return chrDiff;
			}
			return (o2.stopPos - o1.stopPos);
		}
	};
	
	public static GenomeRegion parseRegion (String regionStr) throws AnalysisException {
		// String must be of the form "<chromosome>:<startPos>-<endPos>"
		String[] parts = regionStr.split(":");
		if (parts.length != 2) {
			throw new AnalysisException("Error parsing genomic region '" + regionStr + "': one ':' separator is required.");			
		}
		String chr = parts[0];
		if (chr == null) {
			throw new AnalysisException("Error parsing genomic region '" + regionStr + "': invalid chromosome specified.");			
		}
		String boundaryStr = parts[1];
		String[] boundaryPartStr = boundaryStr.split("-");
		if (boundaryPartStr.length != 2) {
			throw new AnalysisException("Error parsing genomic region '" + regionStr + "': one '-' separator is required.");			
		}
		int startPos, endPos;
		try {
			startPos = Integer.parseInt(boundaryPartStr[0].trim());
			endPos = Integer.parseInt(boundaryPartStr[1].trim());
		} catch (Exception e) {
			throw new AnalysisException("Error parsing start/end positions of genomic region '" + regionStr + "'");			
		}
		if ((startPos <= 0) || (startPos > endPos)) {
			throw new AnalysisException("Error parsing genomic region '" + regionStr + "': invalid boundary positions.");			
		}
		GenomeRegion r = new GenomeRegion(chr, startPos, endPos);
		return r;
	}

}
