package org.cggh.bam.breakpoint;

import java.util.*;

public class Join {
	
	public static final int POLY_T  = 0;
	public static final int POLY_A  = 1;
	public static final int POLY_AT = 2;

	public static final Join[] ALL_JOINS = {
		new Join(POLY_T,  "polyT",  "TTTTTTTTTTTT",     "(T){12,}"),
		new Join(POLY_A,  "polyA",  "AAAAAAAAAAAA",     "(A){12,}"),
		new Join(POLY_AT, "polyAT", "ATATATATATATATAT", "(AT){8,}"),
	};
	
	
	static private HashMap<String,Join> joinTable;
	static {
		joinTable = new HashMap<String,Join>();
		for (int i = 0; i < ALL_JOINS.length; i++) {
			Join j = ALL_JOINS[i];
			joinTable.put(j.id, j);
		}
	}
	
	
	public static Join getJoin (int type) {
		return ALL_JOINS[type];
	}
	
	public static Join getJoin (String id) {
		return joinTable.get(id);
	}
	
	public static Join findJoin (String testString) {
		if (testString == null) return null;
		for (int i = 0; i < ALL_JOINS.length; i++) {
			if (testString.contains(ALL_JOINS[i].match)) {
				return ALL_JOINS[i];
			}
		}
		return null;
	}
	
	static public Join getReverseComplementJoin (int joinType) {
		switch(joinType) {
		case POLY_T:
			return ALL_JOINS[POLY_A];
		case POLY_A:
			return ALL_JOINS[POLY_T];
		case POLY_AT:
			return ALL_JOINS[POLY_AT];
		}
		return null;
	}
	
	/* ***************************************************************************
	 * 
	 */
	private int     type;
	private String  id;
	private String  match;
	private String  regex;

	public Join(int type, String id, String match, String regex) {
		this.type = type;
		this.id = id;
		this.match = match;
		this.regex = regex;
	}
	
    public int getType() {
		return type;
	}

    public String getId() {
		return id;
	}

	public String getMatch() {
		return match;
	}

	public String getRegex() {
		return regex;
	}

}
