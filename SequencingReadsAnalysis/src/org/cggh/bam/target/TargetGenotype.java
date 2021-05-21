package org.cggh.bam.target;

public class TargetGenotype {
	
	protected String ntGenotype;
	protected Target target;
	
	public TargetGenotype() {
		this(null, null);
	}
	
	public TargetGenotype(String ntGenotype, Target target) {
		this.ntGenotype = ntGenotype;
		this.target = target;
	}
	
	public String getNtGenotype () {
		return ntGenotype;
	}

	public Target getTarget () {
		return target;
	}
	
	public boolean isReferenceAllele () {
	    if ((target instanceof AlignmentTarget) && isValidGenotype()) {
	    	String refSeq = ((AlignmentTarget)target).getTargetRefSeq();
		    return (ntGenotype.equals(refSeq));
	    }
    	return false;
	}
	
	public boolean isValidGenotype() {
		return ((ntGenotype != null) && (target != null));
	}

	public static class NoCoverageTargetGenotype extends TargetGenotype {};
	
	public static class LowQualityTargetGenotype extends TargetGenotype {};
}

