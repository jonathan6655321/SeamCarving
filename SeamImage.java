import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;

import javax.imageio.ImageIO;

public class SeamImage {
	private static final double ENTROPY_WEIGHT = -0.01;
	private static final double EDGES_WEIGHT = 6;
	private static final boolean FILL_ENLARGE = true;
	private static final int NUMBER_OF_CORES = Runtime.getRuntime().availableProcessors();
	private static final int INSERTED_SEAM_COLOR = 0x7f6464;
	private BufferedImage Image;
	private double[][] edgeAndEntropyMatrix, entropyMatrix;
	private int[][] edgeMatrix, grayscaleMatrix, grayscale9X9BlurMatrix;
	private int[][] RGBMatrix;
	private boolean changed = false;

	// constructors:
	public SeamImage(String imageFileName) {
		Image = loadImage(imageFileName);
		RGBMatrix = convertImageToRGBMatrix(Image);
		updateEverythingFromRGB();
	}

	public SeamImage(BufferedImage image) {
		Image = image;
		RGBMatrix = convertImageToRGBMatrix(Image);
		updateEverythingFromRGB();
	}

	public SeamImage(SeamImage image) {
		this.Image = image.Image;
		this.RGBMatrix = Matrix.clone(image.RGBMatrix);
		this.edgeMatrix = Matrix.clone(image.edgeMatrix);
		this.entropyMatrix = Matrix.clone(image.entropyMatrix);
		this.grayscaleMatrix = Matrix.clone(image.grayscaleMatrix);
		this.grayscale9X9BlurMatrix = Matrix.clone(image.grayscale9X9BlurMatrix);
		this.edgeAndEntropyMatrix = Matrix.clone(image.edgeAndEntropyMatrix);
	}

	private static int[][] convertImageToRGBMatrix(BufferedImage image) {
		int numberOfRows = image.getHeight();
		int numberOfColumns = image.getWidth();
		int[][] RGBMatrix = new int[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				RGBMatrix[i][j] = image.getRGB(j, i);
			}
		}
		return RGBMatrix;
	}

	private void updateEverythingFromRGB() {
		edgeMatrix = calculateEdgeMatrix(RGBMatrix);
		grayscaleMatrix = calculateGrayscaleMatrix(RGBMatrix);
		grayscale9X9BlurMatrix = calculateGrayscale9x9BlockAvarageMatrix(grayscaleMatrix);
		entropyMatrix = calculateEntropyMatrix(grayscaleMatrix, grayscale9X9BlurMatrix);
		edgeAndEntropyMatrix = calculateEdgeAndEntropyMatrix(entropyMatrix, edgeMatrix);
	}

	private BufferedImage loadImage(String fileName) {
		try {
			return ImageIO.read(new File(fileName));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	// getters:
	public int getHeight() {
		return RGBMatrix.length;
	}

	public int getWidth() {
		return RGBMatrix[0].length;
	}

	public BufferedImage getImage() {
		if (changed) {
			changed = false;
			Image = Matrix.createBufferImageFromRGBMatrix(RGBMatrix);
		}
		return Matrix.createBufferImageFromRGBMatrix(RGBMatrix);
	}

	public double[][] getEdgeAndEntropyMatrix() {
		return Matrix.clone(this.edgeAndEntropyMatrix);
	}

	public int[][] getEdgeMatrix() {
		return Matrix.clone(this.edgeMatrix);
	}

	public int[][] getRGBMatrix() {
		return Matrix.clone(this.RGBMatrix);
	}

	// rotate:
	public void rotate90right() {
		changed = true;
		RGBMatrix = Matrix.rotateArray(RGBMatrix);
		edgeMatrix = Matrix.rotateArray(edgeMatrix);
		entropyMatrix = Matrix.rotateArray(entropyMatrix);
		grayscale9X9BlurMatrix = Matrix.rotateArray(grayscale9X9BlurMatrix);
		grayscaleMatrix = Matrix.rotateArray(grayscaleMatrix);
		edgeAndEntropyMatrix = Matrix.rotateArray(edgeAndEntropyMatrix);
	}

	// matrix Calculator:
	private static int[][] calculateEdgeMatrix(int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		int[][] edgeMatrix = new int[numberOfRows][numberOfColumns];
		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = 0; col < numberOfColumns; col++) {
					edgeMatrix[row][col] = calculateEdgeValue(row, col, RGBMatrix);
				}
			}
		});
		return edgeMatrix;
	}

	private static double[][] calculateEntropyMatrix(int[][] grayscaleMatrix, int[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double[][] entropyMatrix = new double[numberOfRows][numberOfColumns];

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = 0; col < numberOfColumns; col++) {
					entropyMatrix[row][col] = calculateEntropyValue(row, col, grayscaleMatrix, grayscale9X9BlurMatrix);
				}
			}
		});
		return entropyMatrix;
	}

	private static int[][] calculateGrayscale9x9BlockAvarageMatrix(int[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		int[][] grayscale9X9BlurMatrix = new int[numberOfRows][numberOfColumns];
		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = 0; col < numberOfColumns; col++) {
					grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockSumValue(row, col, grayscaleMatrix);
				}
			}
		});
		return grayscale9X9BlurMatrix;
	}

	private static double[][] calculateEdgeAndEntropyMatrix(double[][] entropyMatrix, int[][] edgeMatrix) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;

		double[][] edgeAndEntropyMatrix = new double[numberOfRows][numberOfColumns];

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = 0; col < numberOfColumns; col++) {
					edgeAndEntropyMatrix[row][col] = EDGES_WEIGHT * edgeMatrix[row][col]
							+ ENTROPY_WEIGHT * entropyMatrix[row][col];
				}
			}
		});
		return edgeAndEntropyMatrix;
	}

	private static int[][] calculateGrayscaleMatrix(int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		int[][] GrayscaleMatrix = new int[numberOfRows][numberOfColumns];

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = 0; col < numberOfColumns; col++) {
					GrayscaleMatrix[row][col] = convertRGBToGrayscaleValue(RGBMatrix[row][col]);
				}
			}
		});
		return GrayscaleMatrix;
	}

	// values calculator:
	private static int calculateEdgeValue(int row, int col, int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		int diff = 0;
		int numberOfNeightbors = 0;
		for (int i = Math.max(row - 1, 0); i < Math.min(row + 2, numberOfRows); i++) {
			for (int j = Math.max(col - 1, 0); j < Math.min(col + 2, numberOfColumns); j++) {
				numberOfNeightbors++;
				diff += calculateEnergyDiffBetweenTwoPixels(RGBMatrix[row][col], RGBMatrix[i][j]);
			}
		}
		return diff / (numberOfNeightbors - 1);
	}

	public static int calculateEnergyDiffBetweenTwoPixels(int pixel1, int pixel2) {
		return Math.abs(((pixel1 >> 16) & 0xff) - ((pixel2 >> 16) & 0xff))
				+ Math.abs(((pixel1 >> 8) & 0xff) - ((pixel2 >> 8) & 0xff))
				+ Math.abs(((pixel1) & 0xff) - ((pixel2) & 0xff));
	}

	private static int convertRGBToGrayscaleValue(int pixel) {
		return (((pixel >> 16) & 0xff) + ((pixel >> 8) & 0xff) + (pixel & 0xff));
	}

	private static int calculateGrayscale9x9BlockSumValue(int row, int col, int[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		int numberOfNeightbors = 0;
		int sum = 0;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				sum += grayscaleMatrix[i][j];
			}
		}
		if (numberOfNeightbors != 81) {
			float complement = sum;
			complement /= numberOfNeightbors;
			sum = (int) (81 * complement);
		}
		return sum;
	}

	private static double calculateEntropyValue(int row, int col, int[][] grayscaleMatrix,
			int[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		int numberOfNeightbors = 0;
		double sum = 0;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				// need to add 1 so there wont be any zeros.
				double pValue = (double)((grayscaleMatrix[i][j] + 1)) / ( grayscale9X9BlurMatrix[i][j]);
				sum += pValue * Math.log(pValue);
			}
		}
		if (numberOfNeightbors != 81) {
			double complement = sum;
			complement /= numberOfNeightbors;
			sum = (81 * complement);
		}
		return sum;
	}

	private static int calculateAvarageColorBasedOnNeightbors(int leftPixel, int rightPixel) {
		int currentPixel;
		if (FILL_ENLARGE) {
			currentPixel = ((((leftPixel >> 16) & 0xff) + ((rightPixel >> 16) & 0xff)) / 2) << 16;
			currentPixel += ((((leftPixel >> 8) & 0xff) + ((rightPixel >> 8) & 0xff)) / 2) << 8;
			currentPixel += ((((leftPixel) & 0xff) + ((rightPixel) & 0xff)) / 2);
			return currentPixel;
		} else {
			return INSERTED_SEAM_COLOR;
		}
	}

	// seam functions:
	public void removeVerticalSeam(int[] seamXValues) {
		changed = true;
		updateRGBMatrix(seamXValues);
		updateEdgeMatrix(seamXValues);
		updateGrayscaleMatrix(seamXValues);
		updateGrayscale9x9BlockAvarageMatrix(seamXValues);
		updateEntropyMatrix(seamXValues);
		updateEdgeAndEntropyMatrix(seamXValues);
	}

	public void enlargeImageHorizontallyByK(int[][] kMinSeams) {
		changed = true;
		int k = kMinSeams.length;
		int[][] resMatrix = new int[RGBMatrix.length][RGBMatrix[0].length + k];

		// insert new seams
		for (int i = 0; i < k; i++) {
			for (int row = 0; row < resMatrix.length; row++) {
				resMatrix[row][kMinSeams[i][row]] = -1;
			}
		}

		// insert old matrix
		for (int row = 0; row < RGBMatrix.length; row++) {
			int offset = 0;
			for (int col = 0; col < RGBMatrix[0].length; col++) {
				while (resMatrix[row][col + offset] == -1) {
					offset++;
				}
				resMatrix[row][col + offset] = RGBMatrix[row][col];
			}
		}

		// calculate averages for new seams:
		int leftPixel, currentPixel, rightPixel;

		for (int[] currentRow : resMatrix) {
			leftPixel = rightPixel = currentRow[1];
			currentPixel = currentRow[0];
			if (currentPixel == -1) {
				currentRow[0] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel, rightPixel);
			}
			for (int col = 2; col < currentRow.length; col++) {
				leftPixel = currentPixel;
				currentPixel = rightPixel;
				rightPixel = currentRow[col];
				if (currentPixel == -1) {
					currentRow[col - 1] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel, rightPixel);
				}
			}
			leftPixel = currentPixel;
			currentPixel = rightPixel;
			rightPixel = leftPixel;
			if (currentPixel == -1) {
				currentRow[currentRow.length - 1] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel,
						rightPixel);
			}
		}

		RGBMatrix = resMatrix;
		updateEverythingFromRGB();

	}

	// matrix updates:
	private void updateRGBMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, RGBMatrix);
	}

	private void updateEdgeMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, edgeMatrix);

		int numberOfRows = edgeMatrix.length;// calculateEdgeValue
		int numberOfColumns = edgeMatrix[0].length;
		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = Math.max(seamXValues[row] - 1, 0); col < Math.min(seamXValues[row] + 1,
						numberOfColumns); col++) {
					edgeMatrix[row][col] = calculateEdgeValue(row, col, RGBMatrix);
				}
			}
		});
	}

	private void updateGrayscaleMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, grayscaleMatrix);
	}

	private void updateGrayscale9x9BlockAvarageMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, grayscale9X9BlurMatrix);
		int numberOfRows = grayscale9X9BlurMatrix.length;
		int numberOfColumns = grayscale9X9BlurMatrix[0].length;

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
						numberOfColumns); col++) {
					grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockSumValue(row, col, grayscaleMatrix);
				}
			}
		});
	}

	private void updateEntropyMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, entropyMatrix);
		int numberOfRows = entropyMatrix.length;
		int numberOfColumns = entropyMatrix[0].length;

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
						numberOfColumns); col++) {
					entropyMatrix[row][col] = calculateEntropyValue(row, col, grayscaleMatrix, grayscale9X9BlurMatrix);
				}
			}
		});
	}

	private void updateEdgeAndEntropyMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, edgeAndEntropyMatrix);
		int numberOfRows = edgeAndEntropyMatrix.length;
		int numberOfColumns = edgeAndEntropyMatrix[0].length;

		int rowsHandaledPerIteration = numberOfRows / NUMBER_OF_CORES;
		IntStream.rangeClosed(0, 1 + numberOfRows / rowsHandaledPerIteration).parallel().forEach(row1 -> {
			int max = Math.min(rowsHandaledPerIteration * (row1 + 1), numberOfRows);
			for (int row = row1 * rowsHandaledPerIteration; row < max; row++) {
				for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
						numberOfColumns - 1); col++) {
					edgeAndEntropyMatrix[row][col] = EDGES_WEIGHT * edgeMatrix[row][col]
							+ ENTROPY_WEIGHT * entropyMatrix[row][col];
				}
			}
		});
	}

	// helper functions for debugging:

	public BufferedImage getImagesGrayscale() {
		return Matrix.createBufferImageFromIntMatrix(grayscaleMatrix);
	}

	public BufferedImage getImagesGrayscaleBlur() {
		return Matrix.createBufferImageFromIntMatrix(grayscale9X9BlurMatrix);
	}

	public BufferedImage getImagesEdges() {
		return Matrix.createBufferImageFromIntMatrix(edgeMatrix);
	}

	public BufferedImage getImagesEntropy() {
		return Matrix.createBufferImageFromDoubleMatrix(entropyMatrix);
	}

	public BufferedImage getImageEdgeAndEntropy() {
		return Matrix.createBufferImageFromDoubleMatrix(edgeAndEntropyMatrix);
	}

}
