package org.cggh.bam;

import java.io.File;


public class Sample {
	private String name;
	private File   bamFile;
	private String bamChrMap;
	
	public Sample(String name, File bamFile, String bamChrMap) {
		this.name = name;
		this.bamFile = bamFile;
		this.bamChrMap = bamChrMap;
	}
	
	public String getName() {
		return name;
	}

	public File getBamFile() {
		return bamFile;
	}

	public String getBamChromosomeMap() {
		return bamChrMap;
	}

}
