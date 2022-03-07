package org.cggh.bam.breakpoint;

import org.cggh.bam.breakpoint.Flank.*;
import org.cggh.common.config.*;
import org.cggh.common.exceptions.*;
import java.io.File;
import java.util.ArrayList;

public class BreakpointConfig extends BaseConfig {
	
	public static final String PROP_PREFIX = "breakpoint.";
	
	protected Breakpoint[]     breakpoints;
	protected Breakpoint[][][] bpByJoinType;
	protected Join[]           joins;
	
	public BreakpointConfig (File configFile) throws AnalysisException  {
		super(configFile);
		joins = Join.ALL_JOINS;
		breakpoints = parseBreakpoints ();
		bpByJoinType = organizeBreakpointsByJoinType (breakpoints);
	}
	
	public String getPrintableDisplay() {
	    return this.toString();    		
    }

	/* **************************************************************
	 * Accessors
	 */
	public Breakpoint[] getBreakpoints () {
		return breakpoints;
	}

	public Breakpoint getBreakpoint (String bpId) {
		return findBreakpoint (bpId, breakpoints);
	}
	
	public Breakpoint[] getLbpForJoinType (int joinType) {
		return bpByJoinType[joinType][0];
	}
	
	public Breakpoint[] getRbpForJoinType (int joinType) {
		return bpByJoinType[joinType][1];
	}
	
	private Breakpoint findBreakpoint (String bpId, Breakpoint[] bps) {
		for (int idx = 0; idx < bps.length; idx++) {
			if (bps[idx].id.equals(bpId)) {
				return (bps[idx]);
			}
		}
		return null;
	}
	
	/* **************************************************************
	 * Cofig file parsing
	 */
	private Breakpoint[] parseBreakpoints () throws AnalysisException {
		String[] bpIds = getStringListProperty(PROP_PREFIX+"bps");
		ArrayList<Breakpoint> bpList = new ArrayList<Breakpoint>();
		
		for (int idx = 0; idx < bpIds.length; idx++) {
			String id = bpIds[idx];
			String prefix = PROP_PREFIX+"bp."+id+".";
			String typeStr = getMandatoryProperty(prefix+"type");
			int type;
			if ("RBP".equals(typeStr)) {
				type = Breakpoint.TYPE_RBP;
			} else if ("LBP".equals(typeStr)) {
				type = Breakpoint.TYPE_LBP;
			} else {
				throw new AnalysisException("Invalid type specified for breakpoint '"+id+"': '"+typeStr+"'");
			}
			String joinId = getMandatoryProperty(prefix+"join");
			Join join = Join.getJoin(joinId);
			if (join == null) {
				throw new AnalysisException("Invalid join specified for breakpoint '"+id+"': '"+joinId+"' does not exists");
			}
			LFlank left =  (LFlank)parseFlank(prefix, join, true);
			RFlank right = (RFlank)parseFlank(prefix, join, false);
			
			bpList.add(new Breakpoint(id, type, join, left, right));
		}
		Breakpoint[] bps = bpList.toArray(new Breakpoint[bpList.size()]);
		return bps;
	}
	
	private Flank parseFlank (String prefix, Join join, boolean isLeft) throws AnalysisException {
		String suffix = isLeft ? "left." : "right.";
		String anchor = getProperty(prefix+suffix+"anchor");
		String test = getProperty(prefix+suffix+"test");
		if (isLeft) {
			return new Flank.LFlank(anchor, test, join);
		}
		return new Flank.RFlank(anchor, test, join);
	}
	
	private Breakpoint[][][] organizeBreakpointsByJoinType (Breakpoint[] breakpoints) throws AnalysisException {
		ArrayList<Breakpoint> lbpList = new ArrayList<Breakpoint>();
		ArrayList<Breakpoint> rbpList = new ArrayList<Breakpoint>();
		Breakpoint[][][] bpByJoinType = new Breakpoint[joins.length][2][];
		for (int idx = 0; idx < joins.length; idx++) {
			int joinType = joins[idx].getType();
			for (int bpIdx = 0; bpIdx < breakpoints.length; bpIdx++) {
				Breakpoint bp = breakpoints[bpIdx];
				if (bp.join.getType() == joinType) {
					if (bp.type == Breakpoint.TYPE_LBP) {
						lbpList.add(bp);
					} else {
						rbpList.add(bp);
					}
				}
			}
			bpByJoinType[idx][0] = lbpList.toArray(new Breakpoint[lbpList.size()]);
			bpByJoinType[idx][1] = rbpList.toArray(new Breakpoint[rbpList.size()]);
			lbpList.clear();
			rbpList.clear();
		}
		return bpByJoinType;
	}
}

