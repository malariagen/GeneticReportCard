package org.cggh.common.sequence;

import java.util.*;

public class SequenceSet {

	private Sequence[] sequences;
	private Sequence[] orderedSequences;
	private String     name;

	public SequenceSet (Sequence[] sequences, String name) {
		if ((sequences == null) || (sequences.length == 0)) {
			throw new IllegalArgumentException ("Empty list of sequences");
		}
		this.sequences = sequences;
		this.name = name;
	}

	public Sequence[] getSequences() {
		return sequences;
	}
	
	public String getName() {
		return name;
	}

	public String[] getSequenceIds() {
		String[] ids = new String[sequences.length];
		for (int i = 0; i < sequences.length; i++) {
			ids[i] = sequences[i].getId();
		}
		return ids;
	}
	
	public Sequence getSequence(int index) {
		return sequences[index];
	}
	
	public int getSequenceCount() {
		return sequences.length;
	}
	
	public String[] findSequencesWithPeptide(String searchString) {
		if (searchString == null) {
			return getSequenceIds();
		}
		
		ArrayList<String> v = new ArrayList<String>();
		for (int i = 0; i < sequences.length; i++) {
			String data = sequences[i].getData();
			if (data.contains(searchString)) {
				v.add(sequences[i].getId());
			}
		}
		String[] ids = new String[v.size()];
		ids = v.toArray(ids);
		return ids;
	}

	public String getDefaultOutputFilename () {
		return getName()+".fasta";
	}
		
	public String getSummary() {
		StringBuffer sb = new StringBuffer();
		sb.append("# Sequence Set\t");
		sb.append(getName());
		sb.append("\n# Sequence count\t");
		sb.append(getSequenceCount());
		return sb.toString();
	}

	
	private Sequence[] getOrderedSequences () {
		if (orderedSequences == null) {
			orderedSequences = (Sequence[])sequences.clone();
			Arrays.sort(orderedSequences, new SequenceIdComparator());
		}
		return orderedSequences;
	}
	

	/* ********************************************************************
	 * Sequence Set comparisons
	 */
	/**
	 * Tests whether another alignment contains exactly the same list of sequences.
	 * The test is done by comparing the ids of the sequences.
	 */
	public boolean equals (Object other) {
		if (other instanceof SequenceSet) {
			SequenceSet a2 = (SequenceSet)other;
			Sequence[] s1 =  getOrderedSequences ();
			Sequence[] s2 =  a2.getOrderedSequences ();
			return Arrays.equals(s1, s2);
		}
		return false;
	}
	
	private class SequenceIdComparator implements Comparator<Sequence> {
		public int compare (Sequence s1, Sequence s2) {
			return (s1.getId().compareTo(s2.getId()));
		}
	}
}
