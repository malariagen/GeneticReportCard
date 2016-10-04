package org.cggh.common.normalization;

import java.util.*;

public interface Rankable {
	
	public void setRankScore(double rankScore);
	
	public double getRankScore();

	public Comparator<Rankable> getRankingSorter ();
	
}
