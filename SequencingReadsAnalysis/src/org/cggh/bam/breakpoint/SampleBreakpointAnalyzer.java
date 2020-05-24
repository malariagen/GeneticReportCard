package org.cggh.bam.breakpoint;

import org.cggh.bam.*;
import org.cggh.common.exceptions.*;
import htsjdk.samtools.*;
//import org.apache.commons.logging.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;


public class SampleBreakpointAnalyzer {
	
	//private static Log log = LogFactory.getLog(org.cggh.common.util.ClassUtilities.getCurrentClassName());
	
	private static final int MIN_FLANK_LENGTH = 10;
	
	private BreakpointConfig config;
	private Sample           sample;
	@SuppressWarnings("unused")
	private File             outFolder;
	
	private SamReaderFactory  samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

	/* ==========================================================
	 * Invocation: single sample
	 * ==========================================================
	 */
	public SampleBreakpointAnalyzer (BreakpointConfig config, Sample sample, File outFolder) throws AnalysisException  {
		this.config = config;
		this.sample = sample;
		this.outFolder = outFolder;
	}
	
	public SampleResult analyzeSample () throws AnalysisException, IOException  {
		ArrayList<ReadResult> rrList = new ArrayList<ReadResult>();
		SamReader samReader = samReaderFactory.open(sample.getBamFile());
		SAMRecordIterator it = samReader.iterator();
		while (it.hasNext()) {
			SAMRecord record = it.next();
			// This will ignore all secondary and supplementary alignments, reads not passing vendor filters, as well as any duplicates
			if (record.getFlags() >= 256) {
				continue;
			}
			
			// Does the read contain a join?
			String sequence = record.getReadString();

			// What Join is it?
			Join readJoin = Join.findJoin(sequence);
		    if (readJoin == null) {
		    	continue;
		    }
			
			// Analyze the sequence as-is
		    ReadResult rr = analyzeSequence (record, readJoin);
		    if (rr != null) {
		    	rrList.add(rr);
		    	continue;
		    }
			
			// Analyze the reverse complement sequence 
			Join revJoin = Join.getReverseComplementJoin (readJoin.getType());
			SAMRecordUtil.reverseComplement(record);
		    rr = analyzeSequence (record, revJoin);
		    if (rr != null) {
		    	rrList.add(rr);
		    }
		}				
		ReadResult[] rResults = rrList.toArray(new ReadResult[rrList.size()]);
		return new SampleResult(sample, rResults);
	}
	
	public ReadResult analyzeSequence (SAMRecord record, Join join) throws AnalysisException, IOException  {
		String sequence = record.getReadString();
		
		// We're matching an RxL (left flank of a right-hand breakpoint)
		Breakpoint[] rbps = config.getRbpForJoinType(join.getType());
		for (int bpIdx = 0; bpIdx < rbps.length; bpIdx++) {
			Breakpoint bp = rbps[bpIdx];
			Matcher m = bp.left.anchorRegex.matcher(sequence);
		    if (!m.find()) {
		    	continue;
		    }
		    int startIdx = m.start();
		    int endIdx = m.end();
		    // How many characters are there in the right flank to be tested?
		    int flankLen = sequence.length() - endIdx;
		    if (flankLen < MIN_FLANK_LENGTH) {
		    	continue;
		    }
			
		    // Found a match with a testable right flank!
		    // For now, just return it.
		    String matchSeq = m.group();
		    String flankSeq = new String(sequence.substring(endIdx));
		    String remainderSeq = new String(sequence.substring(0,startIdx));
		    
		    ReadResult rr = new ReadResult(record, bp, matchSeq, flankSeq, remainderSeq);
		    return rr;
		}
	    return null;
	}
	

	/* ************************************************************
	 * Results container classes
	 * ************************************************************
	 */
	public class ReadResult  {
		SAMRecord  record;
		Breakpoint breakpoint;
		String     matchSeq;
		String     flankSeq;
		String     remainderSeq;
		
		public ReadResult(SAMRecord record, Breakpoint breakpoint, String matchSeq, String flankSeq, String remainderSeq) {
			this.record = record;
			this.breakpoint = breakpoint;
			this.matchSeq = matchSeq;
			this.flankSeq = flankSeq;
			this.remainderSeq = remainderSeq;
		}
	}
	
	public class SampleResult  {
		Sample       sample;
		ReadResult[] rResults;
		
		public SampleResult(Sample sample, ReadResult[] rResults) {
			this.sample = sample;
			this.rResults = rResults;
		}
	}
}
