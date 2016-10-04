package org.cggh.common.threading;

import java.util.*;

import org.apache.commons.logging.*;

public class ParallelExecutableManager {
	
	public static final int DEFAULT_THREAD_COUNT = 8;

	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	
	private int threadCount;
	private ArrayList<Runnable> tasks = new ArrayList<Runnable>();
	
	public ParallelExecutableManager () {
		this(DEFAULT_THREAD_COUNT);
	}
		
	public ParallelExecutableManager (int threadCount){
		this.threadCount = threadCount;
	}
	
	public void addTask (Runnable task) {
		tasks.add(task);
	}
	
	public synchronized Runnable getTask () {
		if (tasks.isEmpty()) {
			return null;
		}
		return tasks.remove(0);
	}
		
	private ParallelExecutableThread[] threadPool;
	
	public void startExecution () {
		// Create a thread pool of n threads, and start them
		threadPool = new ParallelExecutableThread[threadCount];
		for (int i = 0; i < threadPool.length; i++) {
			threadPool[i] = new ParallelExecutableThread(this);
			threadPool[i].start();
		}
	}
		
	public void waitForCompletion () {
		// Wait for the threads to finish
		for (int i = 0; i < threadPool.length; i++) {
			try {
				threadPool[i].join();
			} catch (InterruptedException e) {
				log.error("Exception while waiting for parallel task threads to complete: "+e);
			}
		}
	}
		
	public void executeSynchronously () {
		// Create a thread pool of n threads, and start them
		startExecution ();
		// Wait for the threads to finish
		waitForCompletion();
	}

	private class ParallelExecutableThread extends Thread {
		private ParallelExecutableManager mgr;
		public ParallelExecutableThread (ParallelExecutableManager mgr) {
			this.mgr = mgr;
		}
		
		public void run() {
			while (true) {
				Runnable task = mgr.getTask();
				if (task == null) {
					break;
				}
				task.run();
			}
		}
	}
}
