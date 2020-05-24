package org.cggh.bam;

import java.io.File;


public class Sample {
	private String name;
	private File   bamFile;
	
	public Sample(String name, File bamFile) {
		this.name = name;
		this.bamFile = bamFile;
	}
	
	public String getName() {
		return name;
	}

	public File getBamFile() {
		return bamFile;
	}
}
