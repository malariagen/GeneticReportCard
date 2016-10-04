package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.*;


public class LocusResult {

	Locus            locus;
	SampleLocusResult[]    sampleResults;
	TargetAlleleCounters[] alleleSampleCounters;
	
	public LocusResult(Locus locus, Sample[] samples) throws AnalysisException {
		this.locus = locus;
		this.sampleResults = new SampleLocusResult[samples.length];
		
		Anchor[] anchors = locus.getAnchors();
		this.alleleSampleCounters = new TargetAlleleCounters[anchors.length];
		for (int i = 0; i < anchors.length; i++) {
			this.alleleSampleCounters[i] = new TargetAlleleCounters();
		}
	}
	
	public Locus getConfig() {
		return locus;
	}

	public SampleLocusResult[] getSampleResults() {
		return sampleResults;
	}

	public void setSampleResults(SampleLocusResult[] sampleResults) {
		this.sampleResults = sampleResults;
	}

	public int getTargetCount() {
		return locus.getTargets().length;
	}

	public TargetAlleleCounters[] getAlleleSampleCounters() {
		return alleleSampleCounters;
	}

	public void setAlleleSampleCounters(TargetAlleleCounters[] alleleSampleCounters) {
		this.alleleSampleCounters = alleleSampleCounters;
	}
}
