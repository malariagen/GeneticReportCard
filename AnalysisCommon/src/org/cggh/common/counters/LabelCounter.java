package org.cggh.common.counters;

import java.util.Comparator;

public class LabelCounter {
	
	protected String label;
	protected int    count;
	
	public LabelCounter(String label) {
		this.label = new String(label);
		this.count = 0;
	}
	
	public void reset () {
		this.count = 0;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
    }

	public int getCount() {
		return count;
	}

	public void setCount(int count) {
		this.count = count;
	}
	
	public void add(int addCount) {
		this.count += addCount;
	}
	
	public void increment () {
		this.count++;
	}

	public void decrement () {
		this.count--;
	}

	public String getSummary() {
		return appendCountToLabel(label);
	}

	protected String appendCountToLabel(String label) {
		return label+":"+count;
	}

	public static Comparator<LabelCounter> sorterByCountDesc = new Comparator<LabelCounter>() {
		public int compare(LabelCounter o1, LabelCounter o2) {
			return o2.count - o1.count;
		}			
	};
}
