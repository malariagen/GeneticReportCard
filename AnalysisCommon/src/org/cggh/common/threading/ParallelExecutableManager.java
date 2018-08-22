package org.cggh.common.threading;

import java.util.*;
import org.apache.commons.logging.*;

public class ParallelExecutableManager {
	
	public static final int DEFAULT_THREAD_COUNT = 8;

	private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private boolean isComplete;
	private int threadCount;
	private int maxTasksInQueue;
	private ArrayList<Runnable> tasks = new ArrayList<Runnable>();
	
	public ParallelExecutableManager () {
		this(DEFAULT_THREAD_COUNT);
	}
	
	public ParallelExecutableManager (int threadCount){
		this.threadCount = threadCount;
		this.maxTasksInQueue = Integer.MAX_VALUE;
		this.isComplete = false;
	}
	
	public void setComplete () {
		this.isComplete = true;
	}
	
	public int getThreadCount () {
		return threadCount;
	}
	
	public int getMaxTasksInQueue () {
		return maxTasksInQueue;
	}
	
	public void setMaxTasksInQueue (int maxTasksInQueue) {
		this.maxTasksInQueue = maxTasksInQueue;
	}
		
	public synchronized void addTask (Runnable task) {
		while (tasks.size() >= maxTasksInQueue) {
			try {
				wait();
			} catch (InterruptedException e) {}
		}
		tasks.add(task);
	    notifyAll();
	}
	
	public synchronized Runnable getTask () {
		while (tasks.isEmpty()) {
			if (isComplete) {
				return null;
			}
			try {
				wait();
			} catch (InterruptedException e) {}
		}
		Runnable task = tasks.remove(0);
	    notifyAll();
		return task;
	}

	private ParallelExecutableThread[] threadPool;
	
	public void startExecution () {
		// Create a thread pool of n threads, and start them
		threadPool = new ParallelExecutableThread[threadCount];
		for (int i = 0; i < threadPool.length; i++) {
			threadPool[i] = new ParallelExecutableThread(this);
			threadPool[i].start();
			log.info("Started thread "+i);
		}
	}
		
	public void waitForThreadsCompletion () {
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
		waitForThreadsCompletion();
	}

	private class ParallelExecutableThread extends Thread {
		private ParallelExecutableManager mgr;
		public ParallelExecutableThread (ParallelExecutableManager mgr) {
			this.mgr = mgr;
		}
		
		public void run() {
			//System.out.println("Starting thread");
			while (true) {
				//System.out.println("Getting task");
				Runnable task = mgr.getTask();
				//System.out.println("Got task");
				if (task == null) {
					break;
				}
				//System.out.println("Running task");
				task.run();
				//System.out.println("Completed task");
			}
		}
	}
}
