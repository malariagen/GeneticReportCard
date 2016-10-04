package org.cggh.bam.genotyping;

import org.cggh.common.genome.GenomeRegion;

import htsjdk.samtools.*;

public class MappedRead {

	public static final int MIN_PHRED_SCORE = 20;
	
	public static final int MAPPED = 1;
	public static final int REMAPPED = 2;
	public static final int UNMAPPED = 3;

	private String  id;
	private String  sequence;
	private String  quality;
	private int     startPos;
	private String  samString;
	private boolean isReversed;
	private int     mappingStatus;
	private boolean coversTargets;
	private boolean hasLowQualityTarget;
	
	public MappedRead(SAMRecord record, Locus locusConfig, Anchor anchor, int anchorPos, int mappingStatus) {
		this(record, locusConfig, (anchor.getPos().getPos()-anchorPos), mappingStatus);
	}
	
	public MappedRead(SAMRecord record, Locus locus, int readStartPos, int mappingStatus) {
		this.id = record.getReadName();
		this.sequence = record.getReadString();
		this.quality = record.getBaseQualityString();
		this.startPos = readStartPos;
		this.samString = record.getSAMString();
		if (samString.charAt(samString.length()-1) == '\n') { // Trim line terminator if there is one
			samString = samString.substring(0, samString.length()-1);
		}
		this.isReversed = record.getReadNegativeStrandFlag();
		this.mappingStatus = mappingStatus;
		
		// Does it cover any targets?
		coversTargets = false;
		hasLowQualityTarget = false;
		int endPos = readStartPos + sequence.length() - 1;
		Target[] targets = locus.getTargets();
		for (int i = 0; i < targets.length; i++) {
			GenomeRegion targetRegion = targets[i].getTargetRegion();
			int tStartPos = targetRegion.getStartPos();
			int tEndPos = targetRegion.getStopPos();
			if ((startPos <= tStartPos) && (endPos >= tEndPos)) {
				int tStartOffset = tStartPos - startPos;
				int tLen = 1 + tEndPos - tStartPos;
				String tSeq = sequence.substring(tStartOffset, tStartOffset+tLen);
				if (!tSeq.contains("N")) {
					coversTargets = true;
					int minQ = 1000;
					for (int j = 0; j < tLen; j++) {
						int q = getPhredScore(tStartOffset+j);
						if (q < minQ) {
							minQ = q;
						}
					}
					hasLowQualityTarget = hasLowQualityTarget | (minQ < MIN_PHRED_SCORE);
				}
			}
		}
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

	public boolean isCoversTargets() {
		return coversTargets;
	}

	public boolean isHasLowQualityTarget() {
		return hasLowQualityTarget;
	}

	public int getPhredScore (int offset) {
		int phredCode = quality.charAt(offset);
		return (phredCode - 33);
	}
			
	public boolean coversAtLeastOneTargets () {
		return this.coversTargets;
	}
	
	public boolean hasLowQualityTargetSeq () {
		return this.hasLowQualityTarget;
	}
}

