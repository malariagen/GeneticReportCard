package org.cggh.common.util;

public class ClassUtilities {

	public static String getCurrentClassName() {
		
		// The name of the invoking class
		return Thread.currentThread().getStackTrace()[2].getClassName();
	}
	
}
