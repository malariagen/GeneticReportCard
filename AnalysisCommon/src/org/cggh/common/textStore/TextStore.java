package org.cggh.common.textStore;

import java.io.*;

public abstract class TextStore {
	
	public static final String GZIP_EXTENSION = ".gz";
	
	protected File folder;
	protected String filename;
	
	public TextStore (File folder, String filename) {
		this.folder = folder;
		this.filename = filename;
	}
	
	public String getPath () {
		return new File (folder, filename).getAbsolutePath();
	}

	public abstract File getFile ();
}
