package org.cggh.common.textStore;

import java.io.*;
import java.util.zip.GZIPOutputStream;

import org.cggh.common.textStore.OutputTextStore;

public class GzippedOutputTextStore extends OutputTextStore {
	
	private File f = null;
	private OutputStreamWriter writer = null;
	
	public GzippedOutputTextStore (File folder, String baseFilename) {
		super (folder, baseFilename);
		String filename = baseFilename + GZIP_EXTENSION;
		f = new File (folder, filename);
	}
	
	@Override
	public File getFile () {
		return f;
	}

	@Override
	public Writer getWriter(boolean append) throws IOException {
		if (writer == null) {
			if (append == true) {
				throw new UnsupportedOperationException("SANITY CHECK - GzippedOutputTextStore cannot be used for append operations unless a writer is already created. Make code corrections.");
			}
			GZIPOutputStream gzos = new GZIPOutputStream(new FileOutputStream(f));
			writer = new OutputStreamWriter(gzos);
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
