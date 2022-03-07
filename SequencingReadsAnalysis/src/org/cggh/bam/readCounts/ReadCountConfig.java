package org.cggh.bam.readCounts;

import org.cggh.bam.BaseAnalysisConfig;
import org.cggh.common.exceptions.*;
import java.io.*;

public class ReadCountConfig extends BaseAnalysisConfig {
	
	public static final String PROP_PREFIX = "readCount.";
	
	public ReadCountConfig (File configFile) throws AnalysisException  {
		super(configFile, PROP_PREFIX);
	}
    
	public String getPrintableDisplay() {
	    return super.getPrintableDisplay();    		
    }
}
