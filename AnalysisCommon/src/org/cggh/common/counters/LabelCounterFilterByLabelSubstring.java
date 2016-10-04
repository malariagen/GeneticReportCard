package org.cggh.common.counters;

public class LabelCounterFilterByLabelSubstring implements LabelCounterFilter {
	
	private String badString;
	
	public LabelCounterFilterByLabelSubstring (String badString) {
		this.badString = badString;
	}
	
	public boolean isCounterValid (LabelCounter counter) {
		return (!counter.label.contains(badString));
	}
}