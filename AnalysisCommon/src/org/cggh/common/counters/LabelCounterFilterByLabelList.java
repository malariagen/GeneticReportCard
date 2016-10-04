package org.cggh.common.counters;


public class LabelCounterFilterByLabelList implements LabelCounterFilter {
	
	private String[] labelsAllowed;
	
	public LabelCounterFilterByLabelList (String[] labelsAllowed) {
		this.labelsAllowed = labelsAllowed;
	}
	
	public boolean isCounterValid (LabelCounter counter) {
		String l = counter.getLabel();
		for (int i = 0; i < labelsAllowed.length; i++) {
			if (labelsAllowed[i].equals(l)) {
				return true;
			}
		}
		return false;
	}
}

