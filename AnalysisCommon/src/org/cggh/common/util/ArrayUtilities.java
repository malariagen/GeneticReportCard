package org.cggh.common.util;

import java.lang.reflect.Array;
import java.util.*;

public class ArrayUtilities {
	
	public static int[] toIntArray(List<Integer> list)  {
	    int[] ret = new int[list.size()];
	    int i = 0;
	    for (Integer e : list) {
	        ret[i++] = e.intValue();
	    }
	    return ret;
	}

	public static double[] toDoubleArray(List<Double> list)  {
		double[] ret = new double[list.size()];
	    int i = 0;
	    for (Double e : list) {
	    	ret[i++] = e.doubleValue();
	    }
	    return ret;
	}

	public static double[] floatListToDoubleArray(List<Float> list)  {
		double[] ret = new double[list.size()];
	    int i = 0;
	    for (Float e : list) {
	    	ret[i++] = e.doubleValue();
	    }
	    return ret;
	}

	public static float[] toFloatArray(List<Double> list)  {
		float[] ret = new float[list.size()];
	    int i = 0;
	    for (Double e : list) {
	    	ret[i++] = e.floatValue();
	    }
	    return ret;
	}

	public static float[] doubleListToFloatArray(List<Double> list)  {
		float[] ret = new float[list.size()];
	    int i = 0;
	    for (Double e : list) {
	    	ret[i++] = e.floatValue();
	    }
	    return ret;
	}

	public static <T> T[] concatenateArrays(T[] a, T[] b) {
	    int aLen = a.length;
	    int bLen = b.length;
	    @SuppressWarnings("unchecked")
	    T[] c = (T[]) Array.newInstance(a.getClass().getComponentType(), aLen + bLen);
	    System.arraycopy(a, 0, c, 0, aLen);
	    System.arraycopy(b, 0, c, aLen, bLen);
	    return c;
	}
}
