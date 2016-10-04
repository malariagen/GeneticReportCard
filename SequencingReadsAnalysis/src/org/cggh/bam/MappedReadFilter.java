package org.cggh.bam;

import java.util.*;

public abstract class MappedReadFilter {
	
	public abstract boolean isReadValid(MappedRead read);
	
	public ArrayList<MappedRead> filterReads(ArrayList<MappedRead> mappedReadList) {
		Iterator<MappedRead> it = mappedReadList.iterator();
		ArrayList<MappedRead> discardedReadsList = new ArrayList<MappedRead>();
		while(it.hasNext()) {
			MappedRead read = it.next();
			if (!isReadValid(read)) {
				it.remove();
				discardedReadsList.add(read);
			}
		}
		return discardedReadsList;
	}
}
