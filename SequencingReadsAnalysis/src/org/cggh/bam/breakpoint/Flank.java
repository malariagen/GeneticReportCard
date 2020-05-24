package org.cggh.bam.breakpoint;

import java.util.regex.Pattern;

public abstract class Flank {
	
	protected String  anchor;	
   	protected Pattern anchorRegex;
   	protected String  testString;
	
	public Flank(String anchor, String testString, Join join) {
		this.anchor = anchor;
		this.testString = testString;
		this.anchorRegex = buildAnchorRegex(anchor, join);
	}
	
	protected abstract Pattern buildAnchorRegex (String anchor, Join join);
	
	public Pattern getAnchorRegex() {
		return anchorRegex;
	}

	public String getAnchor() {
		return anchor;
	}

	public String getTestString() {
		return testString;
	}

	public static class LFlank extends Flank {
		public LFlank(String anchor, String testString, Join join) {
			super(anchor, testString, join);
			
		}
		protected Pattern buildAnchorRegex (String anchor, Join join) {
			return Pattern.compile(anchor+join.getRegex());
		}
	}
	
	public static class RFlank extends Flank {
		public RFlank(String anchor, String testString, Join join) {
			super(anchor, testString, join);
		}
		protected Pattern buildAnchorRegex (String anchor, Join join) {
			return Pattern.compile(join.getRegex()+anchor);
		}
	}
}
