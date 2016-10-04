package org.cggh.common.sequence.io;

public interface FileTypeProcessor {
	
	/**
	 * Gets a text description of the files that are processed by this processor
	 * 
	 * @return the text description
	 */
	public abstract String getDescription ();
	
	/**
	 * Gets a list of file extensions that are processed by this processor
	 * 
	 * @return an array containing the file extensions (without '.')
	 */
	public abstract String[] getFileExtensions ();

}
