package org.cggh.common.fileIO;

import org.cggh.common.exceptions.AnalysisException;
import org.cggh.common.textStore.InputTextStore;

import java.io.*;
import java.text.NumberFormat;
import java.util.Arrays;

public abstract class SquareMatrix {
	
	protected String[]   individuals;
	
	public SquareMatrix (String[] individuals) {
		this.individuals = individuals;
	}

	public String[] getIndividuals() {
		return individuals;
	}

	/* ****************************************************************************
	 * 
	 */
	public static class SquareMatrixFloat extends SquareMatrix {
		public  float[][] values;
		
		public SquareMatrixFloat (String[] individuals) {
			super(individuals);
			this.values = new float[individuals.length][individuals.length];
		}

		public float[][] getValues() {
			return values;
		}

		public void setValue (int indIdx1, int indIdx2, float value) {
			values[indIdx1][indIdx2]  = value;
			values[indIdx2][indIdx1]  = value;
		}
		
		public void addValue (int indIdx1, int indIdx2, float value) {
			values[indIdx1][indIdx2]  += value;
			values[indIdx2][indIdx1]  += value;
		}
		
		public float getValue (int indIdx1, int indIdx2) {
			return values[indIdx1][indIdx2];
		}
		
		public void writeToFile (File folder, String filename, NumberFormat floatFormat) throws AnalysisException {
			int bufferSize = 32 * 1024 * 1024;
			String[] colNames = new String[individuals.length+1];
			colNames[0] = "Id";
			for (int i = 0; i < individuals.length; i++) {
				colNames[i+1] = individuals[i];
			}
			TableOutput matrixOut = new TableOutput (folder, filename, colNames, bufferSize);
			if (floatFormat != null) {
				matrixOut.setFloatFormat(floatFormat);
			}
			for (int idx1 = 0; idx1 < individuals.length; idx1++) {
				matrixOut.newRow();
				matrixOut.appendValue(individuals[idx1]);
				for (int idx2 = 0; idx2 < individuals.length; idx2++) {
					matrixOut.appendValue(values[idx1][idx2]);
				}
			}
			matrixOut.close();
		}

		public static SquareMatrixFloat readFromFile (File folder, String filename) throws IOException, AnalysisException {
			File inFile = new File (folder, filename);
			TableInput tif = new TableInput (new InputTextStore (inFile));
			String[] fieldNames = tif.getFieldNames();
			String[] indivIds = Arrays.copyOfRange(fieldNames, 1, fieldNames.length);
			SquareMatrixFloat m = new SquareMatrixFloat (indivIds);
			int rowIdx = 0;
			try {
				while (true) {
					String[] fields = tif.getNextValidLine();
					if (fields == null) {
						break;
					}
					String rowName = fields[0];
					if (!rowName.equals(indivIds[rowIdx])) {
						throw new AnalysisException("SquareMatrix parsing exception: row name mismatch on row #"+rowIdx);
					}
					for (int colIdx = 0; colIdx < indivIds.length; colIdx++) {
						float value = Float.parseFloat(fields[colIdx+1]);
						m.values[rowIdx][colIdx] = value;
						m.values[colIdx][rowIdx] = value;
					}
					rowIdx++;
				}
			} finally {
				tif.close();
			}
			return m;
		}
	}

	/* ****************************************************************************
	 * 
	 */
	public static class SquareMatrixDouble extends SquareMatrix {
		public  double[][] values;
		
		public SquareMatrixDouble (String[] individuals) {
			super(individuals);
			this.values = new double[individuals.length][individuals.length];
		}

		public double[][] getValues() {
			return values;
		}

		public void setValue (int indIdx1, int indIdx2, double value) {
			values[indIdx1][indIdx2]  = value;
			values[indIdx2][indIdx1]  = value;
		}
		
		public void addValue (int indIdx1, int indIdx2, double value) {
			values[indIdx1][indIdx2]  += value;
			values[indIdx2][indIdx1]  += value;
		}
		
		public double getValue (int indIdx1, int indIdx2) {
			return values[indIdx1][indIdx2];
		}
		
		public void writeToFile (File folder, String filename, NumberFormat floatFormat) throws AnalysisException {
			int bufferSize = 32 * 1024 * 1024;
			String[] colNames = new String[individuals.length+1];
			colNames[0] = "Id";
			for (int i = 0; i < individuals.length; i++) {
				colNames[i+1] = individuals[i];
			}
			TableOutput matrixOut = new TableOutput (folder, filename, colNames, bufferSize, TableOutput.COMPRESSION_NONE, TableOutput.LINE_NUMBERS_OFF);
			if (floatFormat != null) {
				matrixOut.setFloatFormat(floatFormat);
			}
			for (int idx1 = 0; idx1 < individuals.length; idx1++) {
				matrixOut.newRow();
				matrixOut.appendValue(individuals[idx1]);
				for (int idx2 = 0; idx2 < individuals.length; idx2++) {
					matrixOut.appendValue(values[idx1][idx2]);
				}
			}
			matrixOut.close();
		}

		public static SquareMatrixDouble readFromFile (File folder, String filename) throws IOException, AnalysisException {
			File inFile = new File (folder, filename);
			TableInput tif = new TableInput (new InputTextStore (inFile));
			String[] fieldNames = tif.getFieldNames();
			String[] indivIds = Arrays.copyOfRange(fieldNames, 1, fieldNames.length);
			SquareMatrixDouble m = new SquareMatrixDouble (indivIds);
			int rowIdx = 0;
			try {
				while (true) {
					String[] fields = tif.getNextValidLine();
					if (fields == null) {
						break;
					}
					String rowName = fields[0];
					if (rowName != indivIds[rowIdx]) {
						throw new AnalysisException("SquareMatrix parsing exception: row name mismatch on row #"+rowIdx);
					}
					for (int colIdx = 0; colIdx < indivIds.length; colIdx++) {
						double value = Double.parseDouble(fields[colIdx+1]);
						m.values[rowIdx][colIdx] = value;
						m.values[colIdx][rowIdx] = value;
					}
					rowIdx++;
				}
			} finally {
				tif.close();
			}
			return m;
		}

	}

	/* ****************************************************************************
	 * 
	 */
	public static class SquareMatrixInt extends SquareMatrix {
		public  int[][] values;
		
		public SquareMatrixInt (String[] individuals) {
			super(individuals);
			this.values = new int[individuals.length][individuals.length];
		}
		
		public int[][] getValues() {
			return values;
		}

		public void setValue (int indIdx1, int indIdx2, int value) {
			values[indIdx1][indIdx2]  = value;
			values[indIdx2][indIdx1]  = value;
		}
		
		public void addValue (int indIdx1, int indIdx2, int value) {
			values[indIdx1][indIdx2]  += value;
			values[indIdx2][indIdx1]  += value;
		}

		public void incrementValue (int indIdx1, int indIdx2) {
			values[indIdx1][indIdx2]++;
			values[indIdx2][indIdx1]++;
		}

		public double getValue (int indIdx1, int indIdx2) {
			return values[indIdx1][indIdx2];
		}
		
		public void writeToFile (File folder, String filename) throws AnalysisException {
			int bufferSize = 32 * 1024 * 1024;
			String[] colNames = new String[individuals.length+1];
			colNames[0] = "Id";
			for (int i = 0; i < individuals.length; i++) {
				colNames[i+1] = individuals[i];
			}
			TableOutput matrixOut = new TableOutput (folder, filename, colNames, bufferSize);
			for (int idx1 = 0; idx1 < individuals.length; idx1++) {
				matrixOut.newRow();
				matrixOut.appendValue(individuals[idx1]);
				for (int idx2 = 0; idx2 < individuals.length; idx2++) {
					matrixOut.appendValue(values[idx1][idx2]);
				}
			}
			matrixOut.close();
		}

		public static SquareMatrixInt readFromFile (File folder, String filename) throws IOException, AnalysisException {
			File inFile = new File (folder, filename);
			TableInput tif = new TableInput (new InputTextStore (inFile));
			String[] fieldNames = tif.getFieldNames();
			String[] indivIds = Arrays.copyOfRange(fieldNames, 1, fieldNames.length);
			SquareMatrixInt m = new SquareMatrixInt (indivIds);
			int rowIdx = 0;
			try {
				while (true) {
					String[] fields = tif.getNextValidLine();
					if (fields == null) {
						break;
					}
					String rowName = fields[0];
					if (rowName != indivIds[rowIdx]) {
						throw new AnalysisException("SquareMatrix parsing exception: row name mismatch on row #"+rowIdx);
					}
					for (int colIdx = 0; colIdx < indivIds.length; colIdx++) {
						int value = Integer.parseInt(fields[colIdx+1]);
						m.values[rowIdx][colIdx] = value;
						m.values[colIdx][rowIdx] = value;
					}
					rowIdx++;
				}
			} finally {
				tif.close();
			}
			return m;
		}
	}
}
