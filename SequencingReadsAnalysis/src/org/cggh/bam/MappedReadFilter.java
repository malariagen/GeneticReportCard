package org.cggh.bam;

import java.util.*;

public abstract class MappedReadFilter {
	
	public abstract boolean isReadValid(Read read);
	
	public ArrayList<Read> filterReads(ArrayList<Read> mappedReadList) {
		Iterator<Read> it = mappedReadList.iterator();
		ArrayList<Read> discardedReadsList = new ArrayList<Read>();
		while(it.hasNext()) {
			Read read = it.next();
			if (!isReadValid(read)) {
				it.remove();
				discardedReadsList.add(read);
			}
		}
		return discardedReadsList;
	}
}
