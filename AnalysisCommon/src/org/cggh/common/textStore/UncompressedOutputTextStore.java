package org.cggh.common.textStore;

import java.io.*;

public class UncompressedOutputTextStore extends OutputTextStore {
	
	File f = null;
	FileWriter writer = null;
	
	public UncompressedOutputTextStore (File file) {
		super (file.getParentFile(), file.getName());
		f = file;
	}

	public UncompressedOutputTextStore (File folder, String filename) {
		super (folder, filename);
		f = new File (folder, filename);
	}
	
	@Override
	public File getFile () {
		return f;
	}

	@Override
	public Writer getWriter(boolean append) throws IOException {
		if (writer == null) {
			writer = new FileWriter(f, append);
		} else {
			if (!append) {
				throw new UnsupportedOperationException("SANITY CHECK - Requested writer for OutputTextStore while one is still open. Make code correction.");
			}
		}
		return writer;
	}
	
	@Override
	public void closeWriter() throws IOException {
		if (writer != null) {
			writer.close();
			writer = null;
		}
	}
}
