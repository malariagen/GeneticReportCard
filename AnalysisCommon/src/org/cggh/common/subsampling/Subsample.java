package org.cggh.common.subsampling;

import java.util.*;


public class Subsample {
	
	private boolean[] selected;
	private int selectedCount;

	private Random rnd = new Random();
	
	public Subsample (int sampleSize, int selectedCount) {
		this.selectedCount = selectedCount;
		this.selected = new boolean[sampleSize];
		
		ArrayList<Integer> idxList = getIndexList (sampleSize);
		for (int s = 0; s < selectedCount; s++) {
			int rndIdx = rnd.nextInt(idxList.size());
			int idx = idxList.remove(rndIdx);
			selected[idx] = true;
		}
	}
	
	private ArrayList<Integer> getIndexList (int sampleSize) {
		ArrayList<Integer> idxList = new ArrayList<Integer>(sampleSize);
		for (int i = 0; i < sampleSize; i++) {
			idxList.add(i);
		}
		return idxList;
	}
	
	public boolean isSelected(int i) {
		return selected[i];
	}

	public int getSampleSize() {
		return selected.length;
	}

	public int getSelectedCount() {
		return selectedCount;
	}
}
