import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

import javax.imageio.ImageIO;

public class SeamImage {
	private static final double ENTROPY_WEIGHT = -200;
	private static final double EDGES_WEIGHT = 1;
	private static final boolean FILL_ENLARGE = true;
	private static final int INSERTED_SEAM_COLOR = 0x7f6464;
	private static final double LOG_1_81 = (double) Math.log((double) (1) / (double) (81));
	private BufferedImage Image;
	private double[][] edgeMatrix, entropyMatrix, grayscaleMatrix, grayscale9X9BlurMatrix, edgeAndEntropyMatrix;
	private int[][] RGBMatrix;
	boolean changed = false;

	/*
	 * int RGB = image.getRGB(j, i); int R = (RGB >> 16) & 0xff; int G = (RGB >>
	 * * 8) & 0xff; int B = (RGB) & 0xff;
	 */

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
		ForkJoinPool forkJoinPool = new ForkJoinPool();
		forkJoinPool.submit(() -> edgeMatrix = calculateEdgeMatrix(RGBMatrix)).join();
		forkJoinPool.submit(() -> grayscaleMatrix = calculateGrayscaleMatrix(RGBMatrix)).join();
		forkJoinPool.submit(() -> grayscale9X9BlurMatrix = calculateGrayscale9x9BlockAvarageMatrix(grayscaleMatrix))
				.join();
		forkJoinPool.submit(() -> entropyMatrix = calculateEntropyMatrix(grayscaleMatrix, grayscale9X9BlurMatrix))
				.join();
		forkJoinPool.submit(() -> edgeAndEntropyMatrix = calculateEdgeAndEntropyMatrix(entropyMatrix, edgeMatrix))
				.join();
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
		return Image;
	}

	public double[][] getEdgeAndEntropyMatrix() {
		return this.edgeAndEntropyMatrix;
	}

	public double[][] getEdgeMatrix() {
		return this.edgeMatrix;
	}

	public int[][] getRGBMatrix() {
		return this.RGBMatrix;
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
	private static double[][] calculateEdgeMatrix(int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] edgeMatrix = new double[numberOfRows][numberOfColumns];

		IntStream.rangeClosed(0, numberOfRows - 1).parallel().forEach(row1 -> {
			final int row = row1;
			IntStream.rangeClosed(0, numberOfColumns - 1).parallel().forEach(col -> {
				edgeMatrix[row][col] = calculateEdgeValue(row, col, RGBMatrix);
			});
		});
		return edgeMatrix;
	}

	private static double[][] calculateEntropyMatrix(double[][] grayscaleMatrix, double[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double[][] entropyMatrix = new double[numberOfRows][numberOfColumns];

		IntStream.rangeClosed(0, numberOfRows - 1).parallel().forEach(row1 -> {
			final int row = row1;
			IntStream.rangeClosed(0, numberOfColumns - 1).parallel().forEach(col -> {
				int numberOfNeightbors = 0;
				double sum = 0;
				double pValue;
				for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
					for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
						numberOfNeightbors++;
						// need to add 1 so there wont be any zeros.
						pValue = (grayscaleMatrix[i][j] + 1) * grayscale9X9BlurMatrix[i][j] / 81;
						sum += pValue * Math.log10(pValue);
					}
				}
				entropyMatrix[row][col] = (-sum) / numberOfNeightbors;
			});
		});
		return entropyMatrix;
	}

	private static double[][] calculateGrayscale9x9BlockAvarageMatrix(double[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double[][] grayscale9X9BlurMatrix = new double[numberOfRows][numberOfColumns];
		IntStream.rangeClosed(0, numberOfRows - 1).parallel().forEach(row1 -> {
			final int row = row1;
			IntStream.rangeClosed(0, numberOfColumns - 1).parallel().forEach(col -> {
				grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockAvarageValuePowerMinusOne(row, col,
						grayscaleMatrix);
			});
		});
		return grayscale9X9BlurMatrix;
	}

	private static double[][] calculateEdgeAndEntropyMatrix(double[][] entropyMatrix, double[][] edgeMatrix) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;

		double[][] edgeAndEntropyMatrix = new double[numberOfRows][numberOfColumns];

		IntStream.rangeClosed(0, numberOfRows - 1).parallel().forEach(row1 -> {
			final int row = row1;
			IntStream.rangeClosed(0, numberOfColumns - 1).parallel().forEach(col -> {
				edgeAndEntropyMatrix[row][col] = EDGES_WEIGHT * edgeMatrix[row][col]
						+ ENTROPY_WEIGHT * entropyMatrix[row][col];
			});
		});
		return edgeAndEntropyMatrix;
	}

	private static double[][] calculateGrayscaleMatrix(int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] GrayscaleMatrix = new double[numberOfRows][numberOfColumns];

		IntStream.rangeClosed(0, numberOfRows - 1).parallel().forEach(row1 -> {
			final int row = row1;
			IntStream.rangeClosed(0, numberOfColumns - 1).parallel().forEach(col -> {
				GrayscaleMatrix[row][col] = convertRGBToGrayscaleValue(RGBMatrix[row][col]);
			});
		});
		return GrayscaleMatrix;
	}

	// values calculator:
	private static double calculateEdgeValue(int row, int col, int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double diff = 0;
		int numberOfNeightbors = 0;
		for (int i = Math.max(row - 1, 0); i < Math.min(row + 2, numberOfRows); i++) {
			for (int j = Math.max(col - 1, 0); j < Math.min(col + 2, numberOfColumns); j++) {
				numberOfNeightbors++;
				diff += calculateEnergyDiffBetweenTwoPixels(RGBMatrix[row][col], RGBMatrix[i][j]);
			}
		}
		return diff / (numberOfNeightbors - 1);
	}

	public static double calculateEnergyDiffBetweenTwoPixels(int pixel1, int pixel2) {
		return Math.abs(((pixel1 >> 16) & 0xff) - ((pixel2 >> 16) & 0xff))
				+ Math.abs(((pixel1 >> 8) & 0xff) - ((pixel2 >> 8) & 0xff))
				+ Math.abs(((pixel1) & 0xff) - ((pixel2) & 0xff));
	}

	private static double convertRGBToGrayscaleValue(int pixel) {
		return (double) (((pixel >> 16) & 0xff) + ((pixel >> 8) & 0xff) + (pixel & 0xff)) / 3;
	}

	private static double calculateGrayscale9x9BlockAvarageValuePowerMinusOne(int row, int col,
			double[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double numberOfNeightbors = 0;
		double sum = 0;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				sum += grayscaleMatrix[i][j];
			}
		}
		return numberOfNeightbors / sum;
	}

	private static double calculateEntropyValue(int row, int col, double[][] grayscaleMatrix,
			double[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		int numberOfNeightbors = 0;
		double sum = 0;
		double pValue;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				// need to add 1 so there wont be any zeros.
				pValue = (grayscaleMatrix[i][j] + 1) * grayscale9X9BlurMatrix[i][j] / 81;
				sum += pValue * Math.log10(pValue);
			}
		}
		return (-sum) / numberOfNeightbors;
	}

	private static int calculateAvarageColorBasedOnNeightbors(int leftPixel, int currentPixel, int rightPixel) {
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
		// int[][] kMinSeams = SeamCarve.getKMinSeams(k,
		// getEdgeAndEntropyMatrix(EnergyType.HoG));
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

		for (int[] currentRow : resMatrix) {// TODO:Parralel
			leftPixel = rightPixel = currentRow[1];
			currentPixel = currentRow[0];
			if (currentPixel == -1) {
				currentRow[0] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel, currentPixel,
						rightPixel);
			}
			for (int col = 2; col < currentRow.length; col++) {
				leftPixel = currentPixel;
				currentPixel = rightPixel;
				rightPixel = currentRow[col];
				if (currentPixel == -1) {
					currentRow[col - 1] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel, currentPixel,
							rightPixel);
				}
			}
			leftPixel = currentPixel;
			currentPixel = rightPixel;
			rightPixel = leftPixel;
			if (currentPixel == -1) {
				currentRow[currentRow.length - 1] = currentPixel = calculateAvarageColorBasedOnNeightbors(leftPixel,
						currentPixel, rightPixel);
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
		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 1, 0); col < Math.min(seamXValues[row] + 1,
					numberOfColumns); col++) {
				edgeMatrix[row][col] = calculateEdgeValue(row, col, RGBMatrix);
			}
		}
	}

	private void updateGrayscaleMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, grayscaleMatrix);
	}

	private void updateGrayscale9x9BlockAvarageMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, grayscale9X9BlurMatrix);
		int numberOfRows = grayscale9X9BlurMatrix.length;
		int numberOfColumns = grayscale9X9BlurMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns); col++) {
				grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockAvarageValuePowerMinusOne(row, col,
						grayscaleMatrix);
			}
		}
	}

	private void updateEntropyMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, entropyMatrix);
		int numberOfRows = entropyMatrix.length;
		int numberOfColumns = entropyMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns); col++) {
				entropyMatrix[row][col] = calculateEntropyValue(row, col, grayscaleMatrix, grayscale9X9BlurMatrix);
			}
		}
	}

	private void updateEdgeAndEntropyMatrix(int[] seamXValues) {
		Matrix.removeSeam(seamXValues, edgeAndEntropyMatrix);
		int numberOfRows = edgeAndEntropyMatrix.length;
		int numberOfColumns = edgeAndEntropyMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns - 1); col++) {
				edgeAndEntropyMatrix[row][col] = EDGES_WEIGHT * edgeMatrix[row][col]
						+ ENTROPY_WEIGHT * entropyMatrix[row][col];
			}
		}
	}

	// helper functions for debugging:

	public BufferedImage getImagesGrayscale() {
		return Matrix.createBufferImageFromDoubleMatrix(grayscaleMatrix);
	}

	public BufferedImage getImagesGrayscaleBlur() {
		return Matrix.createBufferImageFromDoubleMatrix(grayscale9X9BlurMatrix);
	}

	public BufferedImage getImagesEdges() {
		return Matrix.createBufferImageFromDoubleMatrix(edgeMatrix);
	}

	public BufferedImage getImageEntropy() {
		return Matrix.createBufferImageFromDoubleMatrix(entropyMatrix);
	}

	public BufferedImage getImageEdgeAndEntropy() {
		return Matrix.createBufferImageFromDoubleMatrix(edgeAndEntropyMatrix);
	}

}
