package org.cggh.bam.genotyping;

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
	public void setName(String name) {
		this.name = name;
	}
	public File getBamFile() {
		return bamFile;
	}
	public void setBamFile(File bamFile) {
		this.bamFile = bamFile;
	}
}
