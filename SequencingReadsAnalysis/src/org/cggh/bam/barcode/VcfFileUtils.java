package org.cggh.bam.barcode;

import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.zip.GZIPInputStream;


public class VcfFileUtils {
	
	static {
		System.setProperty("samjdk.try_use_intel_deflater","false");
	}

	public static InputStream createVcfFileReader (File vcfFile) throws AnalysisException {
		if (!vcfFile.canRead()) {
			throw new AnalysisException("Cannot read VCF file "+vcfFile.getAbsolutePath());
		}
		boolean gzipped  = vcfFile.getName().endsWith(".vcf.gz");

		InputStream fStream;
		try {
			if (gzipped) {
				FileInputStream fis = new MyFileInputStream(vcfFile);
				fStream = new GZIPInputStream(fis, 2*4*1024*1024);
			} else {
				FileInputStream fis = new FileInputStream(vcfFile);
				fStream = new BufferedInputStream(fis, 2*4*1024*1024);
			}
		} catch (IOException e) {
			throw new AnalysisException("Error opening VCF file "+vcfFile.getAbsolutePath()+": "+e);
		}
		return fStream;
	}
	
	/**
	 * Fixes a problem in which available() returns 0, killing the GZIP stream and closing early.
	 * Returning 1 all the time might not be the most rigorous solution, but seems to work.
	 *
	 */
	public static class MyFileInputStream extends FileInputStream {
		public MyFileInputStream (File f) throws IOException {
			super(f);
		}
		
	    public int available() throws IOException {
	        return 1;
	    }
	}
}
