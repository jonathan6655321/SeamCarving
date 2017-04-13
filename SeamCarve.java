import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Queue;

import javax.imageio.ImageIO;

/* TODO:
 * Enlarge bigger then original*2;
 * to get better result - maybe first scale the picture, and then do seam carve
 * make enlargement use updateSeamMatrix too
 */
public class SeamCarve {
	public static final boolean SHOW_IMAGE = true;
	public static final boolean SAVE_IMAGE = true;

	public static void main(String[] args) {
		String imageFilename, outputImageFilePath;
		int width, height;
		EnergyType eType;

		imageFilename = args[0];
		width = Integer.parseInt(args[1]);
		height = Integer.parseInt(args[2]);
		eType = EnergyType.phaseEnergy(args[3]);
		outputImageFilePath = args[4];

		resizePicture(imageFilename, width, height, eType, outputImageFilePath);
	}

	private static void resizePicture(String imageFilename, int width, int height, EnergyType eType,
			String outputImageFilePath) {

		SeamImage seamImage = new SeamImage(imageFilename);
		int originalWidth, originalHeight;
		originalWidth = seamImage.getWidth();
		originalHeight = seamImage.getHeight();

		double[][] minSeamsMatrix;
		int[] seamXValues;

		if (isLegitInput(originalWidth, originalHeight, width, height, eType)) {
			if (width > originalWidth) {
				seamImage.enlargeImageHorizontallyByK(getKMinSeams(width - originalWidth, seamImage, eType));
			}
			if (height > originalHeight) {
				seamImage.rotate90right();
				seamImage.enlargeImageHorizontallyByK(getKMinSeams(height - originalHeight, seamImage, eType));
				seamImage.rotate90right();
				seamImage.rotate90right();
				seamImage.rotate90right();
			}
			if (width < seamImage.getWidth()) {
				minSeamsMatrix = calculateMinSeamsMatrixByEnergyType(seamImage, eType);
				for (int i = 0; i < originalWidth - width; i++) {
					seamXValues = getMinSeam(minSeamsMatrix);
					seamImage.removeVerticalSeam(seamXValues);
					updateMinSeamsMatrixByEnergyType(seamImage, seamXValues, minSeamsMatrix, eType);
				}
			}
			if (height < originalHeight) {
				seamImage.rotate90right();

				minSeamsMatrix = calculateMinSeamsMatrixByEnergyType(seamImage, eType);
				seamXValues = null;
				for (int i = 0; i < originalHeight - height; i++) {
					seamXValues = getMinSeam(minSeamsMatrix);
					seamImage.removeVerticalSeam(seamXValues);
					updateMinSeamsMatrixByEnergyType(seamImage, seamXValues, minSeamsMatrix, eType);

				}
				seamImage.rotate90right();
				seamImage.rotate90right();
				seamImage.rotate90right();
			}
			if (SAVE_IMAGE) {
				saveImageToFile(seamImage, outputImageFilePath);
			}
			if (SHOW_IMAGE) {
				ImageFrame.displayImage(seamImage);
			}
		}
	}

	private static boolean isLegitInput(int originalWidth, int originalHeight, int width, int height,
			EnergyType eType) {
		if (width <= 2) {
			System.err.println("ERROR: Width argument should be bigger then 2");
			return false;
		}
		if (height <= 2) {
			System.err.println("ERROR: Height argument should be bigger then 2");
			return false;
		}
		if (originalWidth <= 2) {
			System.err.println("ERROR: Image width should be bigger then 2");
			return false;
		}
		if (originalHeight <= 2) {
			System.err.println("ERROR: Image height should be bigger then 2");
			return false;
		}
		if (eType == null) {
			System.err.println("ERROR: Energy Type should be 0\1\2");
			System.err.println("		0: Regular energy without entropy term");
			System.err.println("		1: Regular energy with entropy term");
			System.err.println("		2: Forward energy");
			return false;
		}
		return true;
	}

	public static void saveImageToFile(SeamImage im, String outputFilePath) {
		try {
			ImageIO.write(im.getImage(), "JPG", new File(outputFilePath));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static double[][] calculateMinSeamsMatrixByEnergyType(SeamImage image, EnergyType eType) {
		switch (eType) {
		case EnergyWithoutEntropy:
			int[][] matrix = image.getEdgeMatrix();
			int height = matrix.length, width = matrix[0].length;
			double[][] newMatrix = new double[height][width];
			for (int i = 0; i < height; i++)
				for (int j = 0; j < width; j++)
					newMatrix[i][j] = matrix[i][j];
			return calculateMinSeamsMatrix(newMatrix);
		case EnergyWithEntropy:
			return calculateMinSeamsMatrix(image.getEdgeAndEntropyMatrix());
		case EnergyForwarding:
			return calculateMinSeamsMatrixForwarding(image.getEdgeAndEntropyMatrix(), image.getRGBMatrix());
		}
		return null;
	}

	public static double[][] calculateMinSeamsMatrix(double[][] edgeAndEntropyMatrix) {
		double[][] minSeamsMatrix = matrixCopy(edgeAndEntropyMatrix);
		int edgeMatrixWidth = edgeAndEntropyMatrix[0].length;
		int edgeMatrixHeight = edgeAndEntropyMatrix.length;

		for (int row = 1; row < edgeMatrixHeight; row++) {
			for (int col = 0; col < edgeMatrixWidth; col++) {
				double minPath = Double.MAX_VALUE;
				minPath = Math.min(minPath, minSeamsMatrix[row - 1][col]);
				if (col != 0) {
					minPath = Math.min(minPath, minSeamsMatrix[row - 1][col - 1]);
				}
				if (col != edgeMatrixWidth - 1) {
					minPath = Math.min(minPath, minSeamsMatrix[row - 1][col + 1]);
				}
				minSeamsMatrix[row][col] = minPath + edgeAndEntropyMatrix[row][col];
			}
		}
		return minSeamsMatrix;
	}

	public static double[][] calculateMinSeamsMatrixForwarding(double[][] edgeAndEntropyMatrix, int[][] RGBMatrix) {
		double[][] minSeamsMatrix = matrixCopy(edgeAndEntropyMatrix);
		int edgeMatrixWidth = edgeAndEntropyMatrix[0].length;
		int edgeMatrixHeight = edgeAndEntropyMatrix.length;
		for (int row = 1; row < edgeMatrixHeight; row++) {
			for (int col = 0; col < edgeMatrixWidth; col++) {
				minSeamsMatrix[row][col] = calculateMinSeamsMatrixForwardingValue(col, row, edgeAndEntropyMatrix,
						RGBMatrix, minSeamsMatrix);
			}
		}
		return minSeamsMatrix;
	}

	public static double calculateMinSeamsMatrixForwardingValue(int col, int row, double[][] edgeAndEntropyMatrix,
			int[][] RGBMatrix, double[][] minSeamsMatrix) {
		double minPath = Double.MAX_VALUE;
		double costUp = 0, costLeft, costRight;
		if (col != 0) {
			costLeft = Math.abs(((RGBMatrix[row - 1][col] >> 16) & 0xff) - ((RGBMatrix[row][col - 1] >> 16) & 0xff))
					+ Math.abs(((RGBMatrix[row - 1][col] >> 8) & 0xff) - ((RGBMatrix[row][col - 1] >> 8) & 0xff))
					+ Math.abs(((RGBMatrix[row - 1][col]) & 0xff) - ((RGBMatrix[row][col - 1]) & 0xff));

			if (col + 1 != RGBMatrix[0].length) {
				costUp = Math.abs(((RGBMatrix[row][col + 1] >> 16) & 0xff) - ((RGBMatrix[row][col - 1] >> 16) & 0xff))
						+ Math.abs(((RGBMatrix[row][col + 1] >> 8) & 0xff) - ((RGBMatrix[row][col - 1] >> 8) & 0xff))
						+ Math.abs(((RGBMatrix[row][col + 1]) & 0xff) - ((RGBMatrix[row][col - 1]) & 0xff));
				costLeft += Math
						.abs(((RGBMatrix[row][col - 1] >> 16) & 0xff) - ((RGBMatrix[row][col + 1] >> 16) & 0xff))
						+ Math.abs(((RGBMatrix[row][col - 1] >> 8) & 0xff) - ((RGBMatrix[row][col + 1] >> 8) & 0xff))
						+ Math.abs(((RGBMatrix[row][col - 1]) & 0xff) - ((RGBMatrix[row][col + 1]) & 0xff));
				costRight = Math
						.abs(((RGBMatrix[row - 1][col] >> 16) & 0xff) - ((RGBMatrix[row][col + 1] >> 16) & 0xff))
						+ Math.abs(((RGBMatrix[row - 1][col] >> 8) & 0xff) - ((RGBMatrix[row][col + 1] >> 8) & 0xff))
						+ Math.abs(((RGBMatrix[row - 1][col]) & 0xff) - ((RGBMatrix[row][col + 1]) & 0xff))
						+ Math.abs(((RGBMatrix[row][col + 1] >> 16) & 0xff) - ((RGBMatrix[row][col - 1] >> 16) & 0xff))
						+ Math.abs(((RGBMatrix[row][col + 1] >> 8) & 0xff) - ((RGBMatrix[row][col - 1] >> 8) & 0xff))
						+ Math.abs(((RGBMatrix[row][col + 1]) & 0xff) - ((RGBMatrix[row][col - 1]) & 0xff));
				minPath = Math.min(minPath, minSeamsMatrix[row - 1][col + 1] + costRight);
			}
			minPath = Math.min(minPath, minSeamsMatrix[row - 1][col - 1] + costLeft);
		} else if (col + 1 != RGBMatrix[0].length) {
			costRight = Math.abs(((RGBMatrix[row - 1][col] >> 16) & 0xff) - ((RGBMatrix[row][col + 1] >> 16) & 0xff))
					+ Math.abs(((RGBMatrix[row - 1][col] >> 8) & 0xff) - ((RGBMatrix[row][col + 1] >> 8) & 0xff))
					+ Math.abs(((RGBMatrix[row - 1][col]) & 0xff) - ((RGBMatrix[row][col + 1]) & 0xff));
			minPath = Math.min(minPath, minSeamsMatrix[row - 1][col + 1] + costRight);
		}
		minPath = Math.min(minPath, minSeamsMatrix[row - 1][col] + costUp);
		return minPath + edgeAndEntropyMatrix[row][col];
	}

	// update matrix:
	public static void updateMinSeamsMatrixByEnergyType(SeamImage image, int[] seamXValues, double[][] minSeamsMatrix,
			EnergyType eType) {
		switch (eType) {
		case EnergyWithoutEntropy:
			int[][] matrix = image.getEdgeMatrix();
			int height = matrix.length, width = matrix[0].length;
			double[][] newMatrix = new double[height][width];
			for (int i = 0; i < height; i++)
				for (int j = 0; j < width; j++)
					newMatrix[i][j] = matrix[i][j];
			updateMinSeamsMatrix(minSeamsMatrix, seamXValues, newMatrix);
			break;
		case EnergyWithEntropy:
			updateMinSeamsMatrix(minSeamsMatrix, seamXValues, image.getEdgeAndEntropyMatrix());
			break;
		case EnergyForwarding:
			updateMinSeamsMatrixForwarding(minSeamsMatrix, seamXValues, image.getEdgeAndEntropyMatrix(),
					image.getRGBMatrix());
			break;
		}
	}

	public static void updateMinSeamsMatrix(double[][] minSeamsMatrix, int[] seamXValues,
			double[][] edgeAndEntropyMatrix) {
		Matrix.removeSeam(seamXValues, minSeamsMatrix);

		int edgeMatrixWidth = edgeAndEntropyMatrix[0].length;
		int edgeMatrixHeight = edgeAndEntropyMatrix.length;

		boolean[] currentRow = new boolean[edgeMatrixWidth], nextRow = new boolean[edgeMatrixWidth], swapRow;

		double oldValue, newValue;
		int col1;
		for (int row = 1; row < edgeMatrixHeight; row++) {
			swapRow = currentRow;
			currentRow = nextRow;
			nextRow = swapRow;
			Arrays.fill(nextRow, false);

			col1 = seamXValues[row - 1];
			if (col1 < edgeMatrixWidth) {
				currentRow[col1] = true;
			}
			if (col1 > 0) {
				currentRow[col1 - 1] = true;
			}
			if (col1 < edgeMatrixWidth - 1) {
				currentRow[col1 + 1] = true;
			}
			for (int col = 0; col < edgeMatrixWidth; col++) {
				if (currentRow[col]) {

					oldValue = minSeamsMatrix[row][col];
					newValue = minSeamsMatrix[row - 1][col];
					if (col != 0) {
						newValue = Math.min(newValue, minSeamsMatrix[row - 1][col - 1]);
					}
					if (col != edgeMatrixWidth - 1) {
						newValue = Math.min(newValue, minSeamsMatrix[row - 1][col + 1]);
					}
					newValue += edgeAndEntropyMatrix[row][col];

					if (newValue != oldValue) {
						minSeamsMatrix[row][col] = newValue;
						if (row + 1 < edgeMatrixHeight) {
							nextRow[col] = true;
							if (col != 0) {
								nextRow[col - 1] = true;
							}
							if (col != edgeMatrixWidth - 1) {
								nextRow[col + 1] = true;
							}
						}
					}
				}
			}
		}
	}

	public static void updateMinSeamsMatrixForwarding(double[][] minSeamsMatrix, int[] seamXValues,
			double[][] edgeAndEntropyMatrix, int[][] RGBMatrix) {
		Matrix.removeSeam(seamXValues, minSeamsMatrix);

		int edgeMatrixWidth = edgeAndEntropyMatrix[0].length;
		int edgeMatrixHeight = edgeAndEntropyMatrix.length;

		boolean[] currentRow = new boolean[edgeMatrixWidth], nextRow = new boolean[edgeMatrixWidth], swapRow;

		double oldValue, newValue;
		int col1;
		for (int row = 1; row < edgeMatrixHeight; row++) {
			swapRow = currentRow;
			currentRow = nextRow;
			nextRow = swapRow;
			Arrays.fill(nextRow, false);

			col1 = seamXValues[row - 1];
			if (col1 < edgeMatrixWidth) {
				currentRow[col1] = true;
			}
			if (col1 > 0) {
				currentRow[col1 - 1] = true;
			}
			if (col1 < edgeMatrixWidth - 1) {
				currentRow[col1 + 1] = true;
			}
			for (int col = 0; col < edgeMatrixWidth; col++) {
				if (currentRow[col]) {

					oldValue = minSeamsMatrix[row][col];
					newValue = calculateMinSeamsMatrixForwardingValue(col, row, edgeAndEntropyMatrix, RGBMatrix,
							minSeamsMatrix);

					if (newValue != oldValue) {
						minSeamsMatrix[row][col] = newValue;
						if (row + 1 < edgeMatrixHeight) {
							nextRow[col] = true;
							if (col != 0) {
								nextRow[col - 1] = true;
							}
							if (col != edgeMatrixWidth - 1) {
								nextRow[col + 1] = true;
							}
						}
					}
				}
			}
		}

	}

	/*
	 * // // from left to right going down. // private static double
	 * costLeft(int row, int col, int[][] RGBMatrix) { // double cost =
	 * Math.abs(((RGBMatrix[row - 1][col] >> 16) & 0xff) - //
	 * ((RGBMatrix[row][col - 1] >> 16) & 0xff)) // + Math.abs(((RGBMatrix[row -
	 * 1][col] >> 8) & 0xff) - ((RGBMatrix[row][col // - 1] >> 8) & 0xff)) // +
	 * Math.abs(((RGBMatrix[row - 1][col]) & 0xff) - ((RGBMatrix[row][col - //
	 * 1]) & 0xff)); // // if (col + 1 < RGBMatrix[0].length) // cost +=
	 * Math.abs(((RGBMatrix[row][col - 1] >> 16) & 0xff) - //
	 * ((RGBMatrix[row][col + 1] >> 16) & 0xff)) // +
	 * Math.abs(((RGBMatrix[row][col - 1] >> 8) & 0xff) - ((RGBMatrix[row][col
	 * // + 1] >> 8) & 0xff)) // + Math.abs(((RGBMatrix[row][col - 1]) & 0xff) -
	 * ((RGBMatrix[row][col + // 1]) & 0xff)); // return cost; // } // // //
	 * from right to left going down. // private static double costRight(int
	 * row, int col, int[][] RGBMatrix) { // double cost =
	 * Math.abs(((RGBMatrix[row - 1][col] >> 16) & 0xff) - //
	 * ((RGBMatrix[row][col + 1] >> 16) & 0xff)) // + Math.abs(((RGBMatrix[row -
	 * 1][col] >> 8) & 0xff) - ((RGBMatrix[row][col // + 1] >> 8) & 0xff)) // +
	 * Math.abs(((RGBMatrix[row - 1][col]) & 0xff) - ((RGBMatrix[row][col + //
	 * 1]) & 0xff)); // // if (col - 1 >= 0) // cost +=
	 * Math.abs(((RGBMatrix[row][col + 1] >> 16) & 0xff) - //
	 * ((RGBMatrix[row][col - 1] >> 16) & 0xff)) // +
	 * Math.abs(((RGBMatrix[row][col + 1] >> 8) & 0xff) - ((RGBMatrix[row][col
	 * // - 1] >> 8) & 0xff)) // + Math.abs(((RGBMatrix[row][col + 1]) & 0xff) -
	 * ((RGBMatrix[row][col - // 1]) & 0xff)); // return cost; // } // // //
	 * from above // private static double costUp(int row, int col, int[][]
	 * RGBMatrix) { // if (col - 1 >= 0 && col + 1 < RGBMatrix[0].length) { //
	 * double cost = Math.abs(((RGBMatrix[row][col + 1] >> 16) & 0xff) - //
	 * ((RGBMatrix[row][col - 1] >> 16) & 0xff)) // +
	 * Math.abs(((RGBMatrix[row][col + 1] >> 8) & 0xff) - ((RGBMatrix[row][col
	 * // - 1] >> 8) & 0xff)) // + Math.abs(((RGBMatrix[row][col + 1]) & 0xff) -
	 * ((RGBMatrix[row][col - // 1]) & 0xff)); // return cost; // } else //
	 * return 0; // }
	 */

	public static int[] getMinSeam(double[][] minSeamsMatrix) {
		int numCol = minSeamsMatrix[0].length;
		int numRows = minSeamsMatrix.length;
		int[] seam = new int[numRows];
		seam[numRows - 1] = getMinElementsIndex(minSeamsMatrix[numRows - 1]);

		for (int i = numRows - 2; i >= 0; i--) {
			int prevIndex = seam[i + 1];

			double minVal = minSeamsMatrix[i][prevIndex]; // directly above
			int minIndex = prevIndex;

			if (prevIndex != 0 && minVal > minSeamsMatrix[i][prevIndex - 1]) {
				minVal = minSeamsMatrix[i][prevIndex - 1];
				minIndex--;
			}
			if (prevIndex != numCol - 1 && minVal > minSeamsMatrix[i][prevIndex + 1]) {
				minVal = minSeamsMatrix[i][prevIndex + 1];
				minIndex = prevIndex + 1;
			}

			seam[i] = minIndex;
		}
		return seam;
	}

	public static int[][] getKMinSeams(int k, SeamImage seamImage, EnergyType eType) {
		int numRows = seamImage.getEdgeAndEntropyMatrix().length;
		int numCols = seamImage.getEdgeAndEntropyMatrix()[0].length;
		int[][] kMinSeams = new int[k][numRows];
		SeamImage seamImageClone = new SeamImage(seamImage.getImage());

		boolean[][] isInSeam = new boolean[numRows][numCols];

		for (int i = 0; i < k; i++) {
			kMinSeams[i] = getMinSeam(calculateMinSeamsMatrixByEnergyType(seamImageClone, eType));
			seamImageClone.removeVerticalSeam(kMinSeams[i]);
			for (int row = 0; row < numRows; row++) {
				int offset = 0;
				for (int col = 0; col <= kMinSeams[i][row]; col++) {
					while (isInSeam[row][col + offset]) {
						offset++;
					}
				}
				kMinSeams[i][row] += offset;
				isInSeam[row][kMinSeams[i][row]] = true;
			}
		}

		return calculateTrueInsertedIndex(kMinSeams);
	}

	public static int[][] calculateTrueInsertedIndex(int[][] kMinSeams) {
		int numRows = kMinSeams[0].length;
		for (int i = 0; i < numRows; i++) // iterate over rows
		{
			for (int j = 0; j < kMinSeams.length; j++) // iterate over seams
			{
				for (int k = 0; k < kMinSeams.length; k++) // add 1 to index
															// of seams
															// whose index
															// is larger
															// than that of
															// j seam
				{
					if (k != j && kMinSeams[k][i] >= kMinSeams[j][i]) {
						kMinSeams[k][i]++;
					}
				}
			}
		}

		for (int i = 0; i < numRows; i++) // iterate over rows
		{
			for (int j = 0; j < kMinSeams.length; j++) // iterate over seams
			{
				for (int k = j + 1; k < kMinSeams.length; k++) // add 1 to index
																// of seams
																// whose index
																// is larger
																// than that of
																// j seam
				{
					if (kMinSeams[k][i] == kMinSeams[j][i]) {
						System.out.println("error");
					}
				}
			}
		}

		return kMinSeams;
	}

	public static int getMinElementsIndex(double[] arr) {
		double min = arr[0];
		int j = 0;

		for (int i = 0; i < arr.length; i++) {
			if (arr[i] < min) {
				min = arr[i];
				j = i;
			}
		}
		return j;
	}

	public static double[][] matrixCopy(double[][] original) {
		double[][] newMat = new double[original.length][];
		for (int i = 0; i < original.length; i++)
			newMat[i] = original[i].clone();
		return newMat;
	}

}
