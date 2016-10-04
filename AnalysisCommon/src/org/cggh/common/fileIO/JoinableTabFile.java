package org.cggh.common.fileIO;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.GzippedOutputTextStore;
import org.cggh.common.util.FileUtilities;

import java.io.*;
import java.util.*;


public class JoinableTabFile {
	
	protected Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());

	private GzippedOutputTextStore gzippedStore = null;
	private String        baseFilename;
	private File          outFolder;
	
	private ArrayList<File> pieces = new  ArrayList<File>();
	
	public JoinableTabFile (File folder, String filename) {
		this.gzippedStore = new GzippedOutputTextStore(folder, filename);
		this.outFolder    = folder;
		this.baseFilename = filename;
	}
	
	public File getCurrentPiece () {
		if (pieces.isEmpty()) {
			return null;
		}
		int currIdx = pieces.size() - 1;
		return (pieces.get(currIdx));
	}

	public File createPiece () throws AnalysisException {
		int idx = pieces.size() + 1;
		String pieceFilename = baseFilename + ".tmp." + idx;
		File pieceFile = new File (outFolder, pieceFilename);
		FileUtilities.initializeFile(pieceFile);
		pieces.add(pieceFile);
		return pieceFile;
	}
	
	public void joinPieces () throws AnalysisException {
		log.info("Output file Joining - starting "+gzippedStore.getFile().getAbsolutePath());
		FileUtilities.initializeFile(gzippedStore);
		StringBuffer sb = new StringBuffer(64 * 1024 * 1024);  // Big string buffer for efficiency (~256M)

		int pieceCount = pieces.size();
		BufferedReader[] br = new BufferedReader[pieceCount];
		for (int pieceIdx = 0; pieceIdx < pieceCount; pieceIdx++) {
			File pieceFile = pieces.get(pieceIdx);
			try {
				br[pieceIdx] = new BufferedReader(new FileReader(pieceFile));
			} catch (IOException e) {
				throw new AnalysisException("Error while opening file '" + pieceFile.getAbsolutePath() + "' for assembly: " + e);
			}
		}
		
		boolean start = true;
		String[] linePieces = new String[pieceCount];
		int lineNum = 0;
		try {
			while (true) {
				boolean eof = false;
				boolean hasLine = false;
				for (int pieceIdx = 0; pieceIdx < pieceCount; pieceIdx++) {
					linePieces[pieceIdx] = br[pieceIdx].readLine();
					if (linePieces[pieceIdx] == null) {
						eof = true;
						if (hasLine) {
							throw new AnalysisException("Found early EOF in piece '" + pieces.get(pieceIdx).getAbsolutePath() + "'");							
						}
					} else {
						hasLine = true;
						if (eof) {
							throw new AnalysisException("At EOF, found additional lines in piece '" + pieces.get(pieceIdx).getAbsolutePath() + "'");							
						}
					}
				}
				if (eof) {
					break;
				}
				if (start) {
					start = false;
				} else {
					sb.append('\n');
				}
				for (int pieceIdx = 0; pieceIdx < pieceCount; pieceIdx++) {
					sb.append(linePieces[pieceIdx]);
				}
				FileUtilities.commitIfBufferFull(sb, gzippedStore);
				
				if ((lineNum > 0) && ((lineNum % 10000) == 0)) {
					log.info("Output file Joining - processed " + lineNum + " lines");
				}
				lineNum++;
			}
			FileUtilities.commitIfHasContent(sb, gzippedStore);
			gzippedStore.closeWriter();
			log.info("Output file Joining - Completed (" + lineNum + " lines)");

		} catch (Exception e) {
			throw new AnalysisException("Error while processing file assembly at line " + lineNum + ": " + e);
		} finally {
			for (int pieceIdx = 0; pieceIdx < pieceCount; pieceIdx++) {
				File pieceFile = pieces.get(pieceIdx);
				try {
					br[pieceIdx].close();
					pieceFile.delete();
				} catch (Exception e) {
					throw new AnalysisException("Error closing/deleting file '" + pieceFile.getAbsolutePath() + "': " + e);
				}
			}
		}
	}
}
