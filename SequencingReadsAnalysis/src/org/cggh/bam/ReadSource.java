package org.cggh.bam;

import htsjdk.samtools.SAMRecord;

public class ReadSource {
	
	private String  id;
	private String  sequence;
	private String  quality;
	private String  samString;
	private boolean isReversed;

	public ReadSource(SAMRecord record) {
		this.id = record.getReadName();
		this.sequence = record.getReadString();
		this.quality = record.getBaseQualityString();
		this.isReversed = record.getReadNegativeStrandFlag();
		this.samString = record.getSAMString();
		if (samString.charAt(samString.length()-1) == '\n') { 
			samString = samString.substring(0, samString.length()-1);
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

	public String getSamString() {
		return samString;
	}

	public boolean isReversed() {
		return isReversed;
	}

}
