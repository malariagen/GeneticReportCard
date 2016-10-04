package org.cggh.common.textStore;

import java.io.*;

//import com.jcraft.jzlib.GZIPInputStream;
import java.util.zip.*;

/*
 * These classes solves a problem with GZIPInputStream, which seems to close the underlying stream too early.
 * To prevent this, we intercept close() calls until one is explicitly called on this class.
 * See http://www.mail-archive.com/core-libs-dev@openjdk.java.net/msg16571.html
 */
public class GzipInterceptedInputStream  extends GZIPInputStream {
	
	private InterceptedInputStream ccis;
	
	public GzipInterceptedInputStream (File inputFile) throws IOException {
		super(getInterceptedFileInputStream (inputFile));
		ccis = (InterceptedInputStream)in;
	}
	
	@Override
	public void close() throws IOException  {
		super.close();
		ccis.closeUnderlyingStream();
	}
	
	
	private static InterceptedInputStream getInterceptedFileInputStream (File inputFile) throws FileNotFoundException {
		FileInputStream fis = new FileInputStream(inputFile);
		InterceptedInputStream ccis = new InterceptedInputStream (fis);
		return ccis;
	}

	
	public static class InterceptedInputStream extends FilterInputStream {
		public InterceptedInputStream (InputStream in) {
			super(in);
		}
		
		@Override
		public void close() {
			System.out.println("close()");
		}
		
		public void closeUnderlyingStream() throws IOException {
			in.close();
		}
	}
}


