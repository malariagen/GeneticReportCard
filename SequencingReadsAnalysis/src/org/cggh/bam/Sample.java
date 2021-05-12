package org.cggh.bam;

import java.io.File;


public class Sample {
	private String batch;
	private String name;
	private File   bamFile;
	
	public Sample(String batch, String name, File bamFile) {
		this.batch = batch;
		this.name = name;
		this.bamFile = bamFile;
	}
	
	public String getBatch() {
		return batch;
	}

	public String getName() {
		return name;
	}

	public File getBamFile() {
		return bamFile;
	}
}
