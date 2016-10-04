package org.cggh.common.sequence;

public class Sequence {
	
	private String id;
	private String data;
	
	public Sequence (String id, String data) {
		this.id = id;
		this.data = data;
	}
	
	public String getData() {
		return data;
	}
	
	public String getId() {
		return id;
	}
	
	public String getRegion (int startPos, int endPos) {
		int startIdx = startPos - 1;	// 1-based indexing
		int endIdx = endPos;            // start-end pos inclusive
		String regionData = data.substring(startIdx, endIdx);
		return regionData;
	}

}
