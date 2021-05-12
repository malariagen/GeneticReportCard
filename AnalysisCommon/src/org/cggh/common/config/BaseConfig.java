package org.cggh.common.config;

import org.cggh.common.exceptions.AnalysisException;
import java.io.*;
import java.util.*;


public class BaseConfig {
	
	public static final String CONFIG_IMPORT_PROP     = "analysis.config.import";

	protected Properties configProperties;
	
	public BaseConfig (File configFile) throws AnalysisException  {
		this(loadConfigProperties(configFile));
	}

	public BaseConfig (Properties configProperties) {
		this.configProperties = configProperties;
	}

	public Properties getConfigProperties() {
		return configProperties;
	}
	
	public String getPrintableDisplay() {
		return this.toString();
	}
	
	
	/* ==========================================================
	 * Config Property Parsing
	 * ==========================================================
	 */			
	protected String getProperty (String propName) throws AnalysisException {
		String propValue = configProperties.getProperty(propName);
		if (propValue == null) {
			return null;
		}
		propValue = propValue.trim();
		if (propValue.isEmpty()) {
			return null;
		}
		return propValue;
	}
		
	protected String getMandatoryProperty (String propName) throws AnalysisException {
		String propValue = getProperty (propName);
		if (propValue == null) {
			throw new AnalysisException ("Missing value for property '" + propName +"'");
		}
		return propValue;
	}
	
	// *****************************************************************************
	// Boolean Properties
	// *****************************************************************************
	protected boolean getBooleanProperty (String propName) throws AnalysisException {
		String propValueStr = getProperty(propName);
		return parseBooleanProperty (propName, propValueStr);
	}
	protected boolean getBooleanProperty (String propName, boolean defaultValue) throws AnalysisException {
		boolean propValue = defaultValue;
		String propValueStr = getProperty(propName);
		if (propValueStr != null) {
			propValue = parseBooleanProperty (propName, propValueStr);
		}
		return propValue;
	}
	private boolean parseBooleanProperty (String propName, String propValueStr) throws AnalysisException {
		if (propValueStr == null) {
			throw new AnalysisException ("Missing true/false value for property '" + propName +"'");
		}
		if (Boolean.TRUE.toString().equalsIgnoreCase(propValueStr)) {
			return true;
		}
		if (Boolean.FALSE.toString().equalsIgnoreCase(propValueStr)) {
			return false;
		}
		throw new AnalysisException ("Bad true/false value for property '" + propName +"': "+propValueStr);
	}

	// *****************************************************************************
	// Integer Properties
	// *****************************************************************************
	protected int getIntProperty (String propName) throws AnalysisException {
		String propValueStr = getProperty(propName);
		return parseIntProperty (propName, propValueStr);		
	}
	protected int getIntProperty (String propName, int defaultValue) throws AnalysisException {
		int propValue = defaultValue;
		String propValueStr = getProperty(propName);
		if (propValueStr != null) {
			propValue = parseIntProperty (propName, propValueStr);
		}
		return propValue;		
	}
	private int parseIntProperty (String propName, String propValueStr) throws AnalysisException {
		if (propValueStr == null) {
			throw new AnalysisException ("Missing integer value for property '" + propName +"'");
		}
		try {
			int intValue = Integer.parseInt(propValueStr);
			return intValue;
		} catch (Exception e) {
			throw new AnalysisException ("Bad integer value for property '" + propName +"': "+propValueStr);
		}
	}

	// *****************************************************************************
	// Double Properties
	// *****************************************************************************
	protected double getDoubleProperty (String propName) throws AnalysisException {
		String propValueStr = getProperty(propName);
		return parseDoubleProperty (propName, propValueStr);		
	}
	protected double getDoubleProperty (String propName, double defaultValue) throws AnalysisException {
		double propValue = defaultValue;
		String propValueStr = getProperty(propName);
		if (propValueStr != null) {
			propValue = parseDoubleProperty (propName, propValueStr);
		}
		return propValue;			
	}
	private double parseDoubleProperty (String propName, String propValueStr) throws AnalysisException {
		if (propValueStr == null) {
			throw new AnalysisException ("Missing numeric value for property '" + propName +"'");
		}
		try {
			double intValue = Double.parseDouble(propValueStr);
			return intValue;
		} catch (Exception e) {
			throw new AnalysisException ("Bad numeric value for property '" + propName +"': "+propValueStr);
		}
	}
	
	// *****************************************************************************
	// Double Properties
	// *****************************************************************************
	protected File getFileProperty (String propName, boolean checkExists) throws AnalysisException {
		return (getFileProperty (propName, null, checkExists));
	}
	protected File getFileProperty (String propName, File parent, boolean checkExists) throws AnalysisException {
		String fileName = getMandatoryProperty(propName);
		File f = new File (parent, fileName);
		if (checkExists && !f.exists()) {
			throw new AnalysisException("Could not find file '" + f.getAbsolutePath()+"' specified for property '" + propName +"'");			
		}
		return f;
	}
	protected File getFolderProperty (String propName, boolean checkExists) throws AnalysisException {
		File f = getFileProperty (propName, checkExists);
		if (checkExists && !f.isDirectory()) {
			throw new AnalysisException("File '" + f.getAbsolutePath()+"' specified for property '" + propName +"' is not a directory");			
		}
		return f;
	}
	
	// *****************************************************************************
	// String List Properties
	// *****************************************************************************
	protected String[] getStringListProperty (String propName) throws AnalysisException {
		String propValueStr = getProperty(propName);
		if (propValueStr == null) {
			return null;
		}
		String[] propValues = propValueStr.split(",");
		for (int i = 0; i < propValues.length; i++) {
			propValues[i] = propValues[i].trim();
		}
		return propValues;
	}
	protected String[] getStringListProperty (String propName, String[] defaultValues) throws AnalysisException {
		String[] values = getStringListProperty (propName);
		if (values == null) {
			return defaultValues;
		}
		return values;
	}
	
	// *****************************************************************************
	// Integer List Properties
	// *****************************************************************************
	protected int[] getIntegerListProperty (String propName) throws AnalysisException {
		String[] propValues = getStringListProperty (propName);
		if (propValues == null) {
			return null;
		}
		int[] values = new int[propValues.length];
		for (int i = 0; i < propValues.length; i++) {
			try {
				values[i] = Integer.parseInt(propValues[i].trim());
			} catch (Exception e) {
				throw new AnalysisException("Error parsing item #" + (i+1) + " (value='" + propValues[i] + "') in integer list for property '" + propName +"'");
			}
		}
		return values;
	}
	protected int[] getIntegerListProperty (String propName, int[] defaultValues) throws AnalysisException {
		int[] values = getIntegerListProperty (propName);
		if (values == null) {
			return defaultValues;
		}
		return values;
	}
	
	// *****************************************************************************
	// Double List Properties
	// *****************************************************************************
	protected double[] getDoubleListProperty (String propName) throws AnalysisException {
		String[] propValues = getStringListProperty (propName);
		if (propValues == null) {
			return null;
		}
		double[] values = new double[propValues.length];
		for (int i = 0; i < propValues.length; i++) {
			try {
				values[i] = Double.parseDouble(propValues[i].trim());
			} catch (Exception e) {
				throw new AnalysisException("Error parsing item #" + (i+1) + " (value='" + propValues[i] + "') in double list for property '" + propName +"'");
			}
		}
		return values;
	}
	protected double[] getDoubleListProperty (String propName, double[] defaultValues) throws AnalysisException {
		double[] values = getDoubleListProperty (propName);
		if (values == null) {
			return defaultValues;
		}
		return values;
	}
	

	/* ==========================================================
	 * Utilities
	 * ==========================================================
	 */		
	protected static Properties loadConfigProperties(File configFile) throws AnalysisException {
		Properties configProperties = new Properties();
		try {
			Reader reader = new FileReader(configFile);
			configProperties.load(reader);
		} catch (Exception e) {
			throw new AnalysisException("Error parsing configuration file " + configFile.getAbsolutePath() + ": " + e);
		}
		return configProperties;
	}
	
	/* ==========================================================
	 * Config Printing
	 * ==========================================================
	 */		
	/**
	 * Gets a multi-line concatenated display for a bunch of config objects
	 * 
	 * @param configs
	 * @return
	 */
	public static String getPrintableDisplay(BaseConfig[] configs) {
		StringBuffer sb = new StringBuffer();
		boolean first = true;
		for (int i = 0; i < configs.length; i++) {
			BaseConfig config = configs[i];
			if (config != null) {
				if (first) {
					first = false;
				} else {
					sb.append('\n');					
				}
				sb.append(config.getPrintableDisplay());
			}
		}
		return (sb.toString());
	}
	
	protected static String getPrintableString(int[] intList) {
		if (intList == null) {
			return "null";
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < intList.length; i++) {
			if (i > 0) {
				sb.append(',');
			}
			sb.append(intList[i]);
		}
		return sb.toString();
	}
	
	protected static String getPrintableString(double[] dblList) {
		if (dblList == null) {
			return "null";
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < dblList.length; i++) {
			if (i > 0) {
				sb.append(',');
			}
			sb.append(dblList[i]);
		}
		return sb.toString();
	}
	
	protected static String getPrintableString(String[] stringList) {
		if (stringList == null) {
			return "null";
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < stringList.length; i++) {
			if (i > 0) {
				sb.append(',');
			}
			sb.append(stringList[i]);
		}
		return sb.toString();
	}
}
