package org.cggh.bam.genotyping;

public class SampleLocusResult {
	Sample            sample;
	int               readCount;
	int               misalignCount;
	AlignmentTarget[] alignmentTargets;
	
	public SampleLocusResult(Sample sample, int readCount, int misalignCount, AlignmentTarget[] alignmentTargets) {
		this.sample = sample;
		this.readCount = readCount;
		this.misalignCount = misalignCount;
		this.alignmentTargets = alignmentTargets;
	}
}
