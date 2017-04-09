import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.Arrays;

import javax.imageio.ImageIO;

public class SeamImage {
	private static final double ENTROPY_WEIGHT = 0.5;
	private static final double EDGES_WEIGHT = 0.5;
	private BufferedImage originalImage;
	private double[][] edgeMatrix;
	private double[][] entropyMatrix;
	private double[][] grayscale9X9BlurMatrix;
	private double[][] grayscaleMatrix;
	private double[][] edgeAndEntropyMatrix;
	private EnergyType eType = null;
	private int[][][] RGBMatrix;
	boolean changed = false;

	public SeamImage(String imageFileName) {
		originalImage = loadImage(imageFileName);
		RGBMatrix = convertImageToRGB(originalImage);
		edgeMatrix = calculateEdgeMatrix(RGBMatrix, eType);
		grayscaleMatrix = convertRGBToGrayscale(RGBMatrix);
		grayscale9X9BlurMatrix = calculateGrayscale9x9BlockAvarageMatrix(grayscaleMatrix);
		entropyMatrix = calculateEntropyMatrix(grayscaleMatrix, grayscale9X9BlurMatrix);
		edgeAndEntropyMatrix = calculateEdgeAndEntropyMatrix(entropyMatrix, edgeMatrix);
	}

	private BufferedImage loadImage(String fileName) {
		try {
			return ImageIO.read(new File(fileName));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public BufferedImage getOriginalImage() {
		if (changed) {
			changed = false;
			updateBufferImageFromRGB();
		}
		return originalImage;
	}

	private double[][] getEdgeAndEntropyMatrix(EnergyType eType) {
		if (this.eType != eType) {
			edgeMatrix = calculateEdgeMatrix(RGBMatrix, eType);
			this.eType = eType;
			entropyMatrix = calculateEntropyMatrix(grayscaleMatrix, grayscale9X9BlurMatrix);
			edgeAndEntropyMatrix = calculateEdgeAndEntropyMatrix(entropyMatrix, edgeMatrix);
		}
		return edgeAndEntropyMatrix;

	}

	private static double[][] calculateEntropyMatrix(double[][] grayscaleMatrix, double[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double[][] entropyMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				entropyMatrix[i][j] = calculateEntropyValue(i, j, grayscaleMatrix, grayscale9X9BlurMatrix);
			}
		}
		return entropyMatrix;
	}

	private static double[][] calculateGrayscale9x9BlockAvarageMatrix(double[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		double[][] grayscale9X9BlurMatrix = new double[numberOfRows][numberOfColumns];
		for (int row = 0; row < numberOfRows; row++) {
			for (int col = 0; col < numberOfColumns; col++) {
				grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockAvarageValue(row, col, grayscaleMatrix);
			}
		}
		return grayscale9X9BlurMatrix;
	}

	private static double calculateGrayscale9x9BlockAvarageValue(int row, int col, double[][] grayscaleMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;

		int numberOfNeightbors = 0;
		double sum = 0;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				sum += grayscaleMatrix[i][j];
			}
		}
		return sum / numberOfNeightbors;
	}

	private static double[][] calculateEdgeAndEntropyMatrix(double[][] entropyMatrix, double[][] edgeMatrix) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;

		double[][] edgeAndEntropyMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				edgeAndEntropyMatrix[i][j] = EDGES_WEIGHT * edgeMatrix[i][j] + ENTROPY_WEIGHT * entropyMatrix[i][j];
			}
		}
		return edgeAndEntropyMatrix;
	}

	private static double calculateEntropyValue(int row, int col, double[][] grayscaleMatrix,
			double[][] grayscale9X9BlurMatrix) {
		int numberOfRows = grayscaleMatrix.length;
		int numberOfColumns = grayscaleMatrix[0].length;
		double sum = 0;
		int numberOfNeightbors = 0;
		for (int i = Math.max(row - 4, 0); i < Math.min(row + 5, numberOfRows); i++) {
			for (int j = Math.max(col - 4, 0); j < Math.min(col + 5, numberOfColumns); j++) {
				numberOfNeightbors++;
				double pValue = grayscaleMatrix[i][j] / (9 * grayscale9X9BlurMatrix[i][j]);
				sum += pValue * Math.log(pValue);
			}
		}
		return sum / numberOfNeightbors;
	}

	private static double[][] calculateEdgeMatrix(int[][][] RGBMatrix, EnergyType eType) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] edgeMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				edgeMatrix[i][j] = calculateEdgeValue(i, j, RGBMatrix, eType);
			}
		}
		return edgeMatrix;
	}

	private static double calculateEdgeValue(int row, int col, int[][][] RGBMatrix, EnergyType eType) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		if (eType == EnergyType.HoG) {
			int R = RGBMatrix[row][col][0];
			int G = RGBMatrix[row][col][1];
			int B = RGBMatrix[row][col][2];

			double diff = 0;
			int numberOfNeightbors = 0;
			for (int i = Math.max(row - 1, 0); i < Math.min(row + 2, numberOfRows); i++) {
				for (int j = Math.max(col - 1, 0); j < Math.min(col + 2, numberOfColumns); j++) {
					numberOfNeightbors++;
					diff += Math.abs(R - RGBMatrix[i][j][0]) + Math.abs(G - RGBMatrix[i][j][1])
							+ Math.abs(B - RGBMatrix[i][j][2]);
				}
			}
			return diff / (3 * numberOfNeightbors);
		}
		return -1;
	}

	private static int[][][] convertImageToRGB(BufferedImage image) {
		int numberOfRows = image.getHeight();
		int numberOfColumns = image.getWidth();
		int[][][] RGBMatrix = new int[numberOfRows][numberOfColumns][];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				int RGB = image.getRGB(j, i);
				int R = (RGB >> 16) & 0xff;
				int G = (RGB >> 8) & 0xff;
				int B = (RGB) & 0xff;
				RGBMatrix[i][j] = new int[] { R, G, B };
			}
		}
		return RGBMatrix;
	}

	private static double[][] convertRGBToGrayscale(int[][][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] GrayscaleMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				GrayscaleMatrix[i][j] = (double) ((RGBMatrix[i][j][0] + RGBMatrix[i][j][1] + RGBMatrix[i][j][2])) / 3;
			}
		}
		return GrayscaleMatrix;
	}

	public void removeVerticalSeam(int[] seamXValues) {
		changed = true;
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = seamXValues[row]; col < numberOfColumns - 1; col++) {
				RGBMatrix[row][col] = RGBMatrix[row][col + 1];
			}
			RGBMatrix[row] = Arrays.copyOfRange(RGBMatrix[row], 0, numberOfColumns - 1);
		}
	}
	
	// public void removeHorizontalSeam(int[] seamYValues) {
	// int numberOfRows = RGBMatrix.length;
	// int numberOfColumns = RGBMatrix[0].length;
	//
	// for (int col = 0; col < numberOfColumns; col++) {
	// for (int row = seamYValues[col]; row < numberOfRows - 1; row++) {
	// RGBMatrix[row][col] = RGBMatrix[row + 1][col];
	// }
	// RGBMatrix = Arrays.copyOfRange(RGBMatrix, 0, numberOfRows - 1);
	// }
	// }

	private void updateBufferImageFromRGB() {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;

		BufferedImage image = new BufferedImage(numberOfColumns, numberOfRows, BufferedImage.TYPE_INT_RGB);

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = 0; col < numberOfColumns; col++) {
				int rgb = RGBMatrix[row][col][0];
				rgb = (rgb << 8) + RGBMatrix[row][col][1];
				rgb = (rgb << 8) + RGBMatrix[row][col][2];
				image.setRGB(col, row, rgb);
			}
		}
		originalImage = image;
	}
}
