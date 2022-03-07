package org.cggh.bam;

import org.cggh.common.exceptions.*;
import java.util.*;


public interface ReadsRetriever {
	
	public ArrayList<Read>[] retrieveSampleReads (Sample sample) throws AnalysisException;
	
}
