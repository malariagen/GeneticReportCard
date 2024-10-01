package org.cggh.bam.target;

public class TargetGenotype {
	
	protected String ntGenotype;
	protected Target target;
	private String name = null;
	
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
	
	public boolean isValidGenotype() {
		return ((ntGenotype != null) && (target != null));
	}
	
	public String getName() {
		if (name == null) {
			name = target.getName()+":"+ntGenotype;
		}
		return name;
	}
	
	public String toString() {
		return getName();
	}
	
	public static class NoCoverageTargetGenotype extends TargetGenotype {
		public String getName() {
			return "NoCoverage";
		}
	};
	
	public static class LowQualityTargetGenotype extends TargetGenotype {
		public String getName() {
			return "LowQ";
		}		
	};
}

