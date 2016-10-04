package org.cggh.common.counters;

public class LabelCounterFilterByMinCount implements LabelCounterFilter {
	
	private int minCount;
	
	public LabelCounterFilterByMinCount (int minCount) {
		this.minCount = minCount;
	}
	
	public boolean isCounterValid (LabelCounter counter) {
		return (counter.getCount() >= minCount);
	}
}
