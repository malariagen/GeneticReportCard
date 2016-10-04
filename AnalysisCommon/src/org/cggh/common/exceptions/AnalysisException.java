
package org.cggh.common.exceptions;

public class AnalysisException extends Exception {
	
    public AnalysisException () {}

    public AnalysisException (String msg) {
        super (msg);
    }
    
    public AnalysisException (String msg, Exception e) {
            super (msg, e);
    }
    
    public AnalysisException (Exception e) {
            super (e);
    }
}

