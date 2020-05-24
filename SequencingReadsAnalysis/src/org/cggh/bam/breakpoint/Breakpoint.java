package org.cggh.bam.breakpoint;

import org.cggh.bam.breakpoint.Flank.*;

public class Breakpoint {
	
	public static int TYPE_LBP = 1;
	public static int TYPE_RBP = 2;
	
    String id;
    int    type;
    Join   join;
    LFlank  left;
    RFlank  right;
    
	public Breakpoint(String id, int type, Join join, LFlank left, RFlank right) {
		this.id = id;
		this.type = type;
		this.join = join;
		this.left = left;
		this.right = right;
	}
}
