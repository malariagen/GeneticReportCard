package org.cggh.common.threading;

import java.io.*;
import java.util.*;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.fileIO.BufferedTextOutput;

public class PocThreadingTest {

	private static final int MAX_SNP_COUNT_IN_TASK = 100;
	private static final int MAX_TASKS_IN_QUEUE = 4;
	private static final int SNP_COUNT = 850;
	private static final int THREAD_COUNT = 3;
	private static final int COMPUTE_DELAY = 200;
	
	private DataStorage store;
	private ThreadedTaskManager ttm;
	private DistanceProcessTask currentTask = null;
	
	
	public void execute () throws Exception {
		processBeforePass ();
		for (int snpIdx = 0; snpIdx < SNP_COUNT; snpIdx++) {
			processSnp(snpIdx);
		}
		processAfterPass ();
	}
	
	public void processBeforePass () throws AnalysisException {
		store = new DataStorage(new File ("D:\\temp\\pco.txt"));
		ttm = new ThreadedTaskManager(THREAD_COUNT);
		ttm.setMaxTasksInQueue(MAX_TASKS_IN_QUEUE);
		ttm.startExecution();
	}

	
	public void processSnp(int snpIdx) throws AnalysisException {
		System.out.println("Processing "+snpIdx+" of "+SNP_COUNT);
		if (currentTask == null) {
			currentTask = new DistanceProcessTask();
		}
		currentTask.addSnpData(snpIdx);
		if (currentTask.isFull()) {
			System.out.println("Submitting task");
			ttm.addTask(currentTask);
			System.out.println("Submitted task");
			currentTask = null;
		}
	}
	
	public void processAfterPass () throws AnalysisException {
		// Wait till all jobs threaded in parallel complete
		System.out.println("Completing job");
		if (currentTask != null) {
			ttm.addTask(currentTask);
			currentTask = null;
		}
		ttm.setComplete();
		ttm.waitForThreadsCompletion();
		store.close();
	}
	
	
	private class DistanceProcessTask implements Runnable {
		private int maxSnpCount;
		private int snpCount;
		
		private String[] snpData;
		
		
		public DistanceProcessTask () {
			this.maxSnpCount = MAX_SNP_COUNT_IN_TASK;
			snpData = new String[maxSnpCount];
			snpCount = 0;
		}
		
		public boolean isFull() {
			return (maxSnpCount == snpCount);
		}

		public void addSnpData(int snpIdx) {
			snpData[snpCount] = "SNP "+snpIdx;
			snpCount++;
		}

		Random rnd = new Random();
		
		@Override
		public void run() {
			System.out.println("Starting task");
			for (int snpIdx = 0; snpIdx < snpCount; snpIdx++) {
				System.out.println("Computing: "+snpData[snpIdx]);

				// Wait a little, to simulate computation
				try {
					int variation = (COMPUTE_DELAY / 2) - rnd.nextInt(COMPUTE_DELAY);
					Thread.sleep(COMPUTE_DELAY+variation);
				} catch (InterruptedException e) {
					System.err.println("Sleeping error: "+e);
				}
				try {
					store.storeData (snpData[snpIdx]);
				} catch (AnalysisException e) {
					System.err.println("Writing error: "+e);
				} 
				System.out.println("Written: "+snpData[snpIdx]);
			}
		}
	}
	
	private class DataStorage {
		
		private BufferedTextOutput out;
		private StringBuffer       sb;

		public DataStorage(File file) throws AnalysisException  {
			out = new BufferedTextOutput(file, 1024 * 1024);
			sb = out.getStringBuffer();			
		}
		
		public synchronized void storeData (String data) throws AnalysisException  {
			sb.append(data);
			sb.append('\n');
			out.commitIfBufferFull();
		}
		
		public void close () throws AnalysisException {
			out.close();
		}
	}

	/* ==========================================================
	 * Execution
	 * ==========================================================
	 */
	public static void main(String[] args) {
		PocThreadingTest ctrl = null;
		try {
		    ctrl = new PocThreadingTest();
			ctrl.execute();
		} catch (Exception e) {
			System.err.println("Error: " + e);
			return;
		}
	}

}
