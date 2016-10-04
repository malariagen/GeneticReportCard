
package org.cggh.common.util;

import java.util.*;


public class TextUtilities {

	public static String[] mergeStringLists (String[] list1, String[] list2) {
		String[] result = new String[list1.length+list2.length];
		for (int i = 0; i < list1.length; i++) {
			result[i] = list1[i];
		}
		for (int i = 0; i < list2.length; i++) {
			result[i+list1.length] = list2[i];
		}
		return result;
	}

	public static String stringArrayToString (String[] array) {
		return stringArrayToString (array, ",");
	}

	public static String stringArrayToString (String[] array, String delimiter) {
		if ((array == null) || (array.length == 0)){
			return null;
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < array.length; i++) {
			if (i > 0) {
				sb.append(delimiter);
			}
			sb.append(array[i]);
		}
		return sb.toString();
	}

	public static String stringListToString (List<String> list) {
		return stringListToString (list, ",");
	}

	public static String stringListToString (List<String> list, String delimiter) {
		if ((list == null) || list.isEmpty()) {
			return null;
		}
		StringBuffer sb = new StringBuffer();
		boolean start = true;
		for (String s : list) {
			if (start) {
				start = false;
			} else {
				sb.append(delimiter);
			}
			sb.append(s);
		}
		return sb.toString();
	}

	public static String[] stringToStringArray (String value) {
		return stringToStringArray (value, ",", false);
	}

	public static String[] stringToStringArray (String value, boolean useNull) {
		return stringToStringArray (value, ",",  useNull);
	}

	public static String[] stringToStringArray (String value, String delimiter) {
		return stringToStringArray (value, delimiter, false);
	}
	
	public static String[] stringToStringArray (String value, String delimiter, boolean useNull) {
		StringTokenizer st = new StringTokenizer(value, delimiter);
		String[] result = new String[st.countTokens()];
		for (int i = 0; i < result.length; i++) {
			result[i] = st.nextToken().trim();
			if (useNull && (result[i].equals("null"))) {
				result[i] = null;
			}
		}
		return result;
	}
	
	
	
    public static String[] split (String s, int lineLength) {
        int sLength = s.length();
        int numLines = ((sLength - 1) / lineLength) + 1;
        String[] result = new String[numLines];
        for (int i = 0; i < numLines; i++) {
            int start = i * lineLength;
            int finish = start + lineLength;
            finish = (finish > sLength) ? sLength : finish;
			result[i] = s.substring (start, finish);
		}
		return (result);
	}

    public static String normalizeWhitespace (String s) {
		return processWhitespace(s, false);
	}

    public static String removeWhitespace (String s) {
		return processWhitespace(s, true);
	}

    private static String processWhitespace (String s, boolean remove) {
		StringBuffer buf = new StringBuffer();
		boolean isWhite = false;
		boolean atStart = true;
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (Character.isWhitespace(c)) {
				isWhite = true;
			} else {
				if (atStart) {
					atStart = false;
				} else if (isWhite) {
					if (!remove) {
					    buf.append (' ');
					}
				}
				isWhite = false;
                buf.append (c);
			}
		}
		return buf.toString();
	}
    
    

}
