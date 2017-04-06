package org.cggh.common.fileIO;

import org.cggh.common.exceptions.*;
import java.io.*;
import java.util.*;


public class MessageLogFile {

	private ArrayList<Message> messageList;
	private String[]           headers;

	public MessageLogFile(String[] headers) {
		messageList = new ArrayList<Message>();
		this.headers = headers;
	}
	
	public void addMessage (String[] msgFields) {
		messageList.add(new Message(msgFields));
	}
	
	public boolean isEmpty () throws AnalysisException {
		return messageList.isEmpty();
	}
	
	public void clear () throws AnalysisException {
		messageList.clear();
	}
	
	public void saveFile (File outFolder, String filename) throws AnalysisException {
		TableOutput msgOut = new TableOutput(outFolder, filename, headers, 4096);
		for (Message m : messageList) {
			msgOut.newRow();
			for (String f : m.fields) {
				msgOut.appendValue(f);
			}
		}
		msgOut.close();
	}
	
	public static class Message {
		
		protected String[] fields;

		public Message(String[] fields) {
			this.fields = fields;
		}
	}
}
