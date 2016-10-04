package org.cggh.common.normalization;

import java.util.*;


public class Ranker {

	private ArrayList<Rankable> values;
	
	public Ranker () {
		values = new ArrayList<Rankable>();
	}
	
	public void addValue(Rankable value) {
		if (value != null) {
			values.add(value);
		}
	}
	
	public void rank() {
		if (values.isEmpty()) {
			return;
		}
		
		double numValues = (double)values.size();
		Comparator<Rankable> sorter = values.get(0).getRankingSorter();
		Collections.sort(values, sorter);
		
		for (int i = 0; i < numValues; i++) {
			double rankScore = ((double)i) / numValues;
			values.get(i).setRankScore(rankScore);
		}
	}

	public int getValueCount() {
		return values.size();
	}
}
