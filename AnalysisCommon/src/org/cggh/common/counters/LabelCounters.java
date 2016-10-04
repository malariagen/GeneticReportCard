package org.cggh.common.counters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class LabelCounters {
	
	private HashMap<String,LabelCounter> counterTable;

	public LabelCounters() {
		counterTable = new HashMap<String,LabelCounter>();
	}
	
	public LabelCounter createCounter (String label) {
		LabelCounter counter = new LabelCounter(label);
		counterTable.put(label, counter);
		return counter;
	}
	
	public boolean hasCounter (String label) {
		return counterTable.containsKey(label);
	}
	
	public LabelCounter getCounter (String label) {
		return counterTable.get(label);
	}
	
	public boolean isEmpty () {
		return counterTable.isEmpty();
	}
	
	public int getSize () {
		return counterTable.size();
	}

	public int getTotal () {
		int total = 0;
		for (LabelCounter c : counterTable.values()) {
			total += c.getCount();
		}
		return total;
	}
	
	public void setCount (String label, int count) {
		LabelCounter counter = counterTable.get(label);
		if (counter == null) {
			counter = createCounter(label);
		}
		counter.setCount(count);
	}
	
	public void increment (String label) {
		LabelCounter counter = counterTable.get(label);
		if (counter == null) {
			counter = createCounter(label);
		}
		counter.increment();
	}
	
	public void clear() {
		counterTable.clear();
	}
	
	public void filterCounters (LabelCounterFilter filter) {
		Iterator<String> it = counterTable.keySet().iterator(); 
		while (it.hasNext()) {
		    String label = it.next();
			LabelCounter counter = counterTable.get(label);
			if (!filter.isCounterValid(counter)) {
				it.remove();
			}
		}
	}
	
	public LabelCounter[] getSortedCounters () {
		ArrayList<LabelCounter> counterList = new ArrayList<LabelCounter>();
		for (LabelCounter counter : counterTable.values()) {
			counterList.add(counter);
		}
		Collections.sort(counterList, LabelCounter.sorterByCountDesc);
		LabelCounter[] counters = counterList.toArray(new LabelCounter[counterList.size()]);			
		return counters;
	}
	
	public String getSummary () {
		LabelCounter[] counters = getSortedCounters();
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < counters.length; i++) {
			if (i > 0) {
				sb.append(";");
			}
			sb.append(counters[i].getSummary());
		}
		return sb.toString();
	}

}
