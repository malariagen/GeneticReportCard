package org.cggh.bam;

import org.cggh.common.exceptions.*;
import java.util.regex.*;


public class Anchor {
	
	private int            pos;
	private Pattern        regex;
	private String         anchorPattern;
	
	public Anchor (int pos, String anchorPattern) throws AnalysisException {
		this.pos = pos;
		this.anchorPattern = anchorPattern;
		this.regex = Pattern.compile(anchorPattern);
	}

	public int getPos() {
		return pos;
	}

	public Pattern getRegex() {
		return regex;
	}

	public void setRegex(Pattern regex) {
		this.regex = regex;
	}

	public String getAnchorPattern() {
		return anchorPattern;
	}

	public void setAnchorPattern(String anchorPattern) {
		this.anchorPattern = anchorPattern;
	}
}
