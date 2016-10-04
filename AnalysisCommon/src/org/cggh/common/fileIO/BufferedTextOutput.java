package org.cggh.common.fileIO;

import java.io.File;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.*;
import org.cggh.common.util.*;

public class BufferedTextOutput {
	
	private OutputTextStore outputStore;
	private StringBuffer sb;

	public BufferedTextOutput (File outputFile, int bufferSize) throws AnalysisException {
		this (new UncompressedOutputTextStore(outputFile), new StringBuffer(bufferSize));
	}
	
	public BufferedTextOutput (OutputTextStore outputStore, int bufferSize) throws AnalysisException {
		this (outputStore, new StringBuffer(bufferSize));
	}
	
	public BufferedTextOutput (File outputFile, StringBuffer sb) throws AnalysisException {
		this (new UncompressedOutputTextStore(outputFile), sb);
	}
	
	public BufferedTextOutput (OutputTextStore outputStore, StringBuffer sb) throws AnalysisException {
		this.outputStore = outputStore;
		this.sb = sb;
		FileUtilities.initializeFile(outputStore);
	}
	
	public File getOutputFile() {
		return outputStore.getFile();
	}

	public StringBuffer getStringBuffer() {
		return sb;
	}

	public void commitIfBufferFull() throws AnalysisException {
		FileUtilities.commitIfBufferFull (sb, outputStore);	
	}

	public void commitIfHasContent() throws AnalysisException {
		FileUtilities.commitIfHasContent(sb, outputStore);	
	}

	public void close() throws AnalysisException {
		commitIfHasContent();
		try {
			outputStore.closeWriter();
		} catch (Exception e) {
			throw new AnalysisException ("Error closing writer for "+outputStore.getPath()+":" + e);
		}
	}
}
