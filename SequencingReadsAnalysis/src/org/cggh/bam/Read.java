package org.cggh.bam;


public class Read {
	
	public static final int MIN_PHRED_SCORE = 20;
	
	public static final int MAPPED = 1;
	public static final int ANCHORED = 2;
	public static final int UNMAPPED = 3;

	private String  id;
	private String  sequence;
	private String  quality;
	private int     startPos;
	private String  samString;
	private boolean isReversed;
	private int     mappingStatus;
	
	public static Read createMappedRead (ReadSource readSrc, Locus locus, int readStartPos) {
		return new Read (readSrc, locus, readStartPos, MAPPED);
	}
	
	public static Read createAnchoredRead (ReadSource readSrc, Locus locus, Anchor anchor, int anchorPos) {
		return new Read (readSrc, locus, (anchor.getPos()-anchorPos), ANCHORED);
	}
	
	private Read(ReadSource readSrc, Locus locus, int readStartPos, int mappingStatus) {
		this.id = readSrc.getId();
		this.sequence = readSrc.getSequence();
		this.quality = readSrc.getQuality();
		this.startPos = readStartPos;
		this.isReversed = readSrc.isReversed();
		// Trim SAM string line terminator if there is one
		this.samString = readSrc.getSamString();
		if (samString.charAt(samString.length()-1) == '\n') { 
			samString = samString.substring(0, samString.length()-1);
		}
		setMappingStatus (mappingStatus);
	}
	
	public void unmap () {
		setMappingStatus (UNMAPPED);
		this.startPos = -1;
	}
	
	public void setMappingStatus (int mappingStatus) {
		this.mappingStatus = mappingStatus;
	}
	
	public void updateSequence (String sequence, String quality) {
		this.sequence = sequence;
		this.quality = quality;
	}
	
	public String getId() {
		return id;
	}

	public String getSequence() {
		return sequence;
	}


	public String getQuality() {
		return quality;
	}

	public int getStartPos() {
		return startPos;
	}

	public String getSamString() {
		return samString;
	}

	public boolean isReversed() {
		return isReversed;
	}

	public int getMappingStatus() {
		return mappingStatus;
	}
	
	public int getPhredScore (int offset) {
		int phredCode = quality.charAt(offset);
		return (phredCode - 33);
	}
}

