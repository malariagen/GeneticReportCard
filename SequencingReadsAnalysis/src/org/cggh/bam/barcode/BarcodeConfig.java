package org.cggh.bam.barcode;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import java.io.File;


public class BarcodeConfig extends BaseAnalysisConfig {

	private static final String TASK_PROP_PREFIX = "barcode.";
	
	
	public BarcodeConfig(File configFile) throws AnalysisException {
		super(configFile, TASK_PROP_PREFIX);
	}
	
	
	public String getPrintableDisplay() {
	    return super.getPrintableDisplay();    		
    }
}