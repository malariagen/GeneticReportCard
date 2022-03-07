package org.cggh.bam.faceaway;

import org.cggh.common.genome.*;

import htsjdk.samtools.SAMRecord;

public class Faceaway {

	SAMRecord reverseRead;
	SAMRecord forwardRead;
	GenomeRegion region;
	
	public Faceaway(SAMRecord reverseRead, SAMRecord forwardRead) {
		this.reverseRead = reverseRead;
		this.forwardRead = forwardRead;
		String chr = reverseRead.getContig();
		int regStart = reverseRead.getAlignmentEnd();
		int regEnd = forwardRead.getAlignmentStart();
		this.region = new GenomeRegion(chr, regStart, regEnd);
	}
	
	public GenomeRegion getRegion () {
		return region;
	}

	public int getSize () {
		return region.getSize();
	}

	public SAMRecord getReverseRead() {
		return reverseRead;
	}

	public SAMRecord getForwardRead() {
		return forwardRead;
	}
}
