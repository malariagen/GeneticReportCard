package org.cggh.bam.genotyping;

import org.cggh.common.exceptions.*;
import org.cggh.common.genome.*;
import java.util.regex.*;


public class Anchor {
	
	GenomePosition pos;
	Pattern        regex;
	String         anchorPattern;
	
	public Anchor (String posCoords, String anchorPattern) throws AnalysisException {
		this.pos = GenomePosition.parsePosition(posCoords);
		this.anchorPattern = anchorPattern;
		this.regex = Pattern.compile(anchorPattern);
	}

	public GenomePosition getPos() {
		return pos;
	}

	public void setPos(GenomePosition pos) {
		this.pos = pos;
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
