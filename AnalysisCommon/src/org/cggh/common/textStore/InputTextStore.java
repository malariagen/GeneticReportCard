package org.cggh.common.textStore;

import java.io.*;
import java.util.zip.GZIPInputStream;

import org.cggh.common.exceptions.AnalysisException;

//import org.cggh.pf.data.vcf.VcfFileProcessorFactory.MyFileInputStream;
//import org.cggh.pf.util.TabFileReader;

public class InputTextStore extends TextStore {
	
	public InputTextStore (File file) throws AnalysisException {
		this (file.getParentFile(), file.getName());
	}
	
	public InputTextStore (File folder, String baseFilename) throws AnalysisException {
		super (folder, baseFilename);
		File f = getFile();
		if (f == null) {
			throw new AnalysisException("File " + baseFilename + " not found in folder " + folder.getAbsolutePath());
		}
		if (!f.exists()) {
			throw new AnalysisException("File " + f.getAbsolutePath() + " does not exist.");
		} else if (!f.canRead()) {
			throw new AnalysisException("File " + f.getAbsolutePath() + " cannot be read.");
		}
	}
	
	@Override
	public File getFile () {
		// First try to see if there is an uncompressed version of the file
		File f = new File (folder, filename);
		if (f.exists()) {
			return f;
		}

		// No uncompressed version of the file, try the gzipped version
		f = new File (folder, filename + GZIP_EXTENSION);
		if (f.exists()) {
			return f;
		}
		return null;
	}
	
	/**
	 * Get a Reader for the file; if a GZIPped version is found, the reader handles GZIP compression, otherwise use a plain file reader.
	 * 
	 * @return the Reader o
	 * @throws IOException
	 */
	public Reader getReader() throws IOException {
		
		// First try to see if there is a gzipped viersion of the file
		String gzFilename = filename.endsWith(GZIP_EXTENSION) ? filename : (filename + GZIP_EXTENSION);
		File gzfile = new File (folder, gzFilename);
		if (gzfile.exists()) {
			// We have a GZIP file, open a reader with 1M buffer for efficiency
			//GzipInterceptedInputStream gzis = new GzipInterceptedInputStream (f);
			FileInputStream fis = new MyFileInputStream(gzfile);
			GZIPInputStream gzStream = new GZIPInputStream(fis, 1024*1024);
			Reader fReader = new InputStreamReader(gzStream);
			return fReader;
		}

		// No gzipped version of the file, try the unzipped version
		File f = new File (folder, filename);
		if (f.exists()) {
			// We have a plain file, open a plain text file reader
			Reader fReader = new FileReader(f);
			return fReader;
		}
		
		// No version was found!
		throw new java.io.FileNotFoundException("File " + f.getAbsolutePath() + " was not found, either gzipped or unzipped");
	}
	
	
	/**
	 * Fixes a problem in which available() returns 0, killing the GZIP stream and closing early.
	 * Returning 1 all the time might not be the most rigorous solution, but seems to work.
	 *
	 */
	public static class MyFileInputStream extends FileInputStream {
		public MyFileInputStream (File f) throws FileNotFoundException {
			super(f);
		}
		
	    public int available() throws IOException {
	        return 1;
	    }
	}
}
