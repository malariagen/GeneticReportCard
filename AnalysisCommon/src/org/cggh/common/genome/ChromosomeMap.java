package org.cggh.common.genome;

import org.cggh.common.exceptions.*;
import org.cggh.common.fileIO.*;
import org.cggh.common.textStore.*;
import java.io.*;
import java.util.*;


public class ChromosomeMap {
	
	private static HashMap<String, Integer>       chrIndexTable;
	private static HashMap<String, ChromosomeMap> chrMapTable;
	
	private static ChromosomeMap defaultInstance;
	
	
	public static void initialize (File mapFile) throws AnalysisException {
		
		ColumnFileReader cfr = new ColumnFileReader(new InputTextStore(mapFile));
		String[] colnames = cfr.getColumnNames();
		@SuppressWarnings("unchecked")
		ArrayList<String>[] chrNameLists = new ArrayList[colnames.length];
		for (int i = 0; i < colnames.length; i++) {
			chrNameLists[i] = new ArrayList<String>();
		}
		while (cfr.nextRecord()) {
			String[] values = cfr.getValues();
			for (int i = 0; i < colnames.length; i++) {
				chrNameLists[i].add(values[i].trim());
			}
		}
		cfr.close();
		
		int chrCount = chrNameLists[0].size();
		String[][] chrNames = new String[colnames.length][chrCount];
		for (int i = 0; i < colnames.length; i++) {
			for (int j = 0; j < chrCount; j++) {
				String v = chrNameLists[i].get(j);
				chrNames[i][j] = "-".equals(v) ? null : v;
			}
		}
		
		chrIndexTable = new HashMap<String, Integer>();
		for (int j = 0; j < chrCount; j++) {
			chrIndexTable.put(chrNameLists[0].get(j), j);
		}
		
		chrMapTable = new HashMap<String, ChromosomeMap>();
		defaultInstance = new ChromosomeMap ("default", chrNames[0]);
		chrMapTable.put("default", defaultInstance);
		for (int i = 1; i < colnames.length; i++) {
			ChromosomeMap map = new ChromosomeMap (colnames[i], chrNames[i]);
			chrMapTable.put(colnames[i], map);
		}
	}

	public static String getMappedChromosomeName (int chrIdx, String mapName) {
		ChromosomeMap map = getInstance (mapName);
		if (map == null) {
			return null;
		}
		return map.getMappedChromosomeName (chrIdx);
	}

	public static String getMappedChromosomeName (String chrName, String mapName) {
		ChromosomeMap map = getInstance (mapName);
		if (map == null) {
			return null;
		}
		return map.getMappedChromosomeName (chrName);
	}

	public static ChromosomeMap getDefaultInstance () {
		return defaultInstance;
	}
	
	public static ChromosomeMap getInstance (String mapName) {
		if ((mapName == null) || (mapName.equals("default"))) {
			return defaultInstance;
		}
		return chrMapTable.get(mapName);
	}

	
	@SuppressWarnings("unused")
	private String   mapName;
	private String[] chrNames;
	
	private ChromosomeMap (String mapName, String[] chrNames) {
		this.mapName = mapName;
		this.chrNames = chrNames;
	}
	
	public String[] getAllChromosomeNames() {
		return chrNames;
	}

	public String getMappedChromosomeName (int chrIdx) {
		return chrNames[chrIdx];
	}

	public String getMappedChromosomeName (String chrName) {
		Integer chrIdxObj = chrIndexTable.get(chrName);
		if (chrIdxObj == null) {
			return null;
		}
		return getMappedChromosomeName (chrIdxObj.intValue());
	}
}
