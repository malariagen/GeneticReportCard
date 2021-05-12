package org.cggh.common.util;

import java.io.*;
import java.util.ArrayList;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.OutputTextStore;

public class FileUtilities {
	
	public static String getFilenameRoot (File f) {
		String fullname = f.getName();
		String ext = getFileExtension (f);
		if (ext == null) {
			return fullname;
		}
		if (ext.length() == fullname.length()) {
			return "";
		}
		int nameLength = fullname.length() - ext.length() - 1;
		return fullname.substring(0, nameLength);
	}

	public static String getFileExtension (File f) {
		String fullname = f.getName();
		int dotPos = fullname.lastIndexOf('.');
		if (dotPos < 0) {
			return null;
		}
		return fullname.substring(dotPos+1);
	}
	
    public static File checkFolder(File parentFolder, String folderName, boolean createIfMissing) throws AnalysisException {
		File folder = new File (parentFolder, folderName);
		if (createIfMissing && !folder.exists()) {
			folder.mkdir();
		}
		if (!folder.isDirectory()) {
			throw new AnalysisException(folder.getAbsolutePath() + " is not a valid folder.");
		}
		return folder;
	}
	
    public static File checkFolder(File folder, boolean createIfMissing) throws AnalysisException {
		File parentFolder = folder.getParentFile();
		String folderName = parentFolder.getName();
		return checkFolder(parentFolder, folderName, createIfMissing);
	}
	
    public static File getHashedSubfolder (File rootFolder, String fileName) throws AnalysisException {
        int hashcode = fileName.hashCode();
        int firstDirCode = hashcode & 255;
        String firstDirName = String.format("%03d", firstDirCode);
        File firstDir = checkFolder(rootFolder, firstDirName, true);
        int secondDirCode = (hashcode >> 8) & 255;
        String secondDirName = String.format("%03d", secondDirCode);
        File secondDir = checkFolder(firstDir, secondDirName, true);
        return secondDir;
    }
    
    /* ********************************************************************************************
     * PLAIN file version of file output routines
     */
    public static void initializeFile (File file) throws AnalysisException {    	
		FileWriter writer;
		try {
			writer = new FileWriter(file);
			writer.close();
		} catch (IOException e) {
		    throw new AnalysisException ("Error initializing empty file '" + file.getAbsolutePath() + "': "+e);				
		}
    }

    public static void commitIfBufferFull (StringBuffer sb, File currFile) throws AnalysisException {
		int fullMark = (int)(0.98 * (double)sb.capacity());
		boolean isBufferFull = (sb.length() >= fullMark); // Avoid new allocation
		if (isBufferFull) {
			commitIfHasContent (sb, currFile);
		}
	}
	
    public static void commitIfHasContent (StringBuffer sb, File currFile) throws AnalysisException {
		if (sb.length() > 0) {
			try {
				FileUtilities.appendFileContent(sb.toString(), currFile);
			} catch (IOException e) {
			    throw new AnalysisException ("Error committing data to file '" + currFile.getAbsolutePath() + "': "+e);				
			} finally {
				sb.setLength(0);					
			}
		}
	}

    public static void appendFileContent (String content, File file) throws IOException {
    	writeContentToFile (content, file, true);
    }

    public static void writeFileContent (String content, File file) throws IOException {
    	writeContentToFile (content, file, false);
    }

    private static void writeContentToFile (String content, File file, boolean append) throws IOException {
		if (content == null) {
			return;
		}
		FileWriter writer = new FileWriter(file, append);
		writer.write(content);
		writer.close();
    }
    
    
    /* ********************************************************************************************
     * TextStore version of file output routines (can be gzipped)
     */
    public static void initializeFile (OutputTextStore otStore) throws AnalysisException {
    	try {
			otStore.getWriter();   // Just open the file
		} catch (IOException e) {
			throw new AnalysisException("Error initializing OutputTextStore: "+e);
		}
    }
    
    public static void commitIfBufferFull (StringBuffer sb, OutputTextStore otStore) throws AnalysisException {
		int fullMark = (int)(0.98 * (double)sb.capacity());
		boolean isBufferFull = (sb.length() >= fullMark); // Avoid new allocation
		if (isBufferFull) {
			commitIfHasContent (sb, otStore);
		}
	}
	
    public static void commitIfHasContent (StringBuffer sb, OutputTextStore otStore) throws AnalysisException {
		if (sb.length() > 0) {
			try {
				String content = sb.toString();
				Writer writer = otStore.getWriter(true);
				writer.write(content);
			} catch (IOException e) {
			    throw new AnalysisException ("Error committing data to file '" + otStore.getFile().getAbsolutePath() + "': "+e);				
			} finally {
				sb.setLength(0);					
			}
		}
	}
    
    
    /* ********************************************************************************************
     * Read list of entries from file (handy with SNP lists, etc.). We ignore (skip) blank lines by default.
     */
    public static String[] readStringListFromFile (File inputFile) throws AnalysisException {
    	return readStringListFromFile (inputFile, false);
    }
    
    public static String[] readStringListFromFile (File file, boolean errorOnBlankEntry) throws AnalysisException {
    	
    	ArrayList<String> entryList = new ArrayList<String>();
		BufferedReader br = null;
		int lineNum = 0;
		try {
			br = new BufferedReader(new FileReader(file));
			String line;
			while ((line = br.readLine()) != null) {
				lineNum++;
				line = line.trim();
				
				if (line.isEmpty()) {
					if (errorOnBlankEntry) {
						throw new AnalysisException("Found empty line in SNP list file" + file.getAbsolutePath() + " at line " + lineNum);				
					} else {
						continue;
					}
				} else {
					entryList.add(line);
				}
			}
		} catch (IOException e) {
 			throw new AnalysisException("Error reading SNP list file" + file.getAbsolutePath() + " at line " + lineNum + ": "+e);				
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				throw new AnalysisException("Error closing SNP list file" + file.getAbsolutePath() + ": "+e);				
			}
		}
		
		String[] entries = entryList.toArray(new String[entryList.size()]);
		return entries;
	}
    
    public static File[] readFileListFromFile (File inputFile) throws AnalysisException {
        return readFileListFromFile (inputFile, false);
    }
   
    public static File[] readFileListFromFile (File inputFile, boolean errorOnBlankEntry) throws AnalysisException {
    	String[] filenames = readStringListFromFile (inputFile, errorOnBlankEntry);
    	File[] files = new File[filenames.length];
    	for (int i = 0; i < filenames.length; i++) {
    		files[i] = new File(filenames[i]);
    	}
    	return files;
    }
}
