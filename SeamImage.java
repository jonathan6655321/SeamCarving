import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class SeamImage {
	private BufferedImage originalImage;
	private double[][] edgeMatrix;
	private EnergyType eType;
	private int[][][] RGBMatrix;
	private double[][] GrayscaleMatrix;

	public SeamImage(String imageFileName) {
		originalImage = loadImage(imageFileName);
		RGBMatrix = convertImageToRGB(originalImage);
		GrayscaleMatrix = convertRGBToGrayscale(RGBMatrix);
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
		return originalImage;
	}
	public double[][] getEdgeMatrix(EnergyType eType) {
		edgeMatrix = calculateEdgeMatrix(RGBMatrix, eType);
		this.eType = eType;
		return edgeMatrix;
	}

	private static double[][] addEntropy(double[][] GrayscaleMatrix, double[][] edgeMatrix) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;

		double[][] entropyMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 3; i < numberOfRows - 4; i++) {
			for (int j = 3; j < numberOfColumns - 4; j++) {
				//entropyMatrix[i][j] = calculateEdgeValue(i, j, RGBMatrix, eType);// TODO
																					// calc
			}
		}
		addEdgesToMatrix(entropyMatrix, 4);
		return entropyMatrix;
	}

	private static void applyEntropyToEdges(double[][] entropyMatrix, double[][] edgeMatrix) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				edgeMatrix[i][j] -= entropyMatrix[i][j];
			}
		}
	}


	private static double calculateEntropyValue(int row, int col, double[][] GrayscaleMatrix) {
		double diff = 0;
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				//diff += Math.abs(R - RGBMatrix[row + i][col + j][0]) + Math.abs(G - RGBMatrix[row + i][col + j][1])
					//	+ Math.abs(B - RGBMatrix[row + i][col + j][2]);
			}
		}
		return diff / 8;
	}

	private static double[][] calculateEdgeMatrix(int[][][] RGBMatrix, EnergyType eType) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] edgeMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 1; i < numberOfRows - 1; i++) {
			for (int j = 1; j < numberOfColumns - 1; j++) {
				edgeMatrix[i][j] = calculateEdgeValue(i, j, RGBMatrix, eType);
			}
		}
		addEdgesToMatrix(edgeMatrix, 1);
		return edgeMatrix;
	}

	private static void addEdgesToMatrix(double[][] edgeMatrix, int numberOfPixelToAdd) {
		int numberOfRows = edgeMatrix.length;
		int numberOfColumns = edgeMatrix[0].length;
		if (numberOfColumns > numberOfPixelToAdd) {
			for (int k = 0; k < numberOfPixelToAdd; k++) {
				for (int i = k; i < numberOfRows - k; i++) {
					edgeMatrix[i][k] = edgeMatrix[i][k + 1];
					edgeMatrix[i][numberOfColumns - 1 - k] = edgeMatrix[i][numberOfColumns - 2 - k];
				}
			}
		}

		if (numberOfRows > numberOfPixelToAdd) {
			for (int k = 0; k < numberOfPixelToAdd; k++) {
				for (int j = k; j < numberOfColumns - k; j++) {
					edgeMatrix[k][j] = edgeMatrix[k + 1][j];
					edgeMatrix[numberOfRows - 1 - k][j] = edgeMatrix[numberOfRows - 2 - k][j];
				}
			}
		}
	}

	private static double calculateEdgeValue(int row, int col, int[][][] RGBMatrix, EnergyType eType) {
		if (eType == EnergyType.HoG) {
			int R = RGBMatrix[row][col][0];
			int G = RGBMatrix[row][col][1];
			int B = RGBMatrix[row][col][2];

			double diff = 0;
			int numberOfNeightbors = 0;
			for (int i = -1; i < 2; i++) {
				for (int j = -1; j < 2; j++) {
					numberOfNeightbors++;
					diff += Math.abs(R - RGBMatrix[row + i][col + j][0]) + Math.abs(G - RGBMatrix[row + i][col + j][1])
							+ Math.abs(B - RGBMatrix[row + i][col + j][2]);
				}
			}
			return diff / 24;
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

}
