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
		RGBMatrix = convertImageToRGBMatrix(originalImage);
		updateEverythingFromRGB();
	}
	
	private void updateEverythingFromRGB(){
		edgeMatrix = calculateEdgeMatrix(RGBMatrix, eType);
		grayscaleMatrix = convertRGBToGrayscaleMatrix(RGBMatrix);
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

	public double[][] getEdgeAndEntropyMatrix(EnergyType eType) {
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
		return (-sum) / numberOfNeightbors;
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

	private static int[][][] convertImageToRGBMatrix(BufferedImage image) {
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

	private static double[][] convertRGBToGrayscaleMatrix(int[][][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		double[][] GrayscaleMatrix = new double[numberOfRows][numberOfColumns];

		for (int i = 0; i < numberOfRows; i++) {
			for (int j = 0; j < numberOfColumns; j++) {
				GrayscaleMatrix[i][j] = convertRGBToGrayscaleValue(i, j, RGBMatrix);
			}
		}
		return GrayscaleMatrix;
	}

	private static double convertRGBToGrayscaleValue(int row, int col, int[][][] RGBMatrix) {
		return (double) ((RGBMatrix[row][col][0] + RGBMatrix[row][col][1] + RGBMatrix[row][col][2])) / 3;
	}

	public void removeVerticalSeam(int[] seamXValues) {
		changed = true;
		updateRGBMatrix(seamXValues);
		updateEdgeMatrix(seamXValues);
		updateGrayscaleMatrix(seamXValues);
		updateGrayscale9x9BlockAvarageMatrix(seamXValues);
		updateEntropyMatrix(seamXValues);
		updateEdgeAndEntropyMatrix(seamXValues);
	}

	private void updateRGBMatrix(int[] seamXValues) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = seamXValues[row]; col < numberOfColumns - 1; col++) {
				RGBMatrix[row][col] = RGBMatrix[row][col + 1];
			}
		}
		removeLastColomn(RGBMatrix);
	}

	private void updateEdgeMatrix(int[] seamXValues) {
		int numberOfRows = edgeMatrix.length;// calculateEdgeValue
		int numberOfColumns = edgeMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 1, 0); col < Math.min(seamXValues[row] + 1,
					numberOfColumns - 1); col++) {
				edgeMatrix[row][col] = calculateEdgeValue(row, col, RGBMatrix, eType);
			}
			for (int col = seamXValues[row] + 2; col < numberOfColumns - 1; col++) {
				edgeMatrix[row][col] = edgeMatrix[row][col + 1];
			}
		}
		removeLastColomn(edgeMatrix);
	}

	private void updateGrayscaleMatrix(int[] seamXValues) {
		int numberOfRows = grayscaleMatrix.length;// calculateEdgeValue
		int numberOfColumns = grayscaleMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 1, 0); col < Math.min(seamXValues[row] + 1,
					numberOfColumns - 1); col++) {
				grayscaleMatrix[row][col] = convertRGBToGrayscaleValue(row, col, RGBMatrix);
			}
			for (int col = seamXValues[row] + 2; col < numberOfColumns - 1; col++) {
				grayscaleMatrix[row][col] = grayscaleMatrix[row][col + 1];
			}
		}
		removeLastColomn(grayscaleMatrix);
	}

	private void updateGrayscale9x9BlockAvarageMatrix(int[] seamXValues) {
		int numberOfRows = grayscale9X9BlurMatrix.length;
		int numberOfColumns = grayscale9X9BlurMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns - 1); col++) {
				grayscale9X9BlurMatrix[row][col] = calculateGrayscale9x9BlockAvarageValue(row, col, grayscaleMatrix);
			}
			for (int col = seamXValues[row] + 6; col < numberOfColumns - 1; col++) {
				grayscale9X9BlurMatrix[row][col] = grayscale9X9BlurMatrix[row][col + 1];
			}
		}
		removeLastColomn(grayscale9X9BlurMatrix);
	}

	private void updateEntropyMatrix(int[] seamXValues) {
		int numberOfRows = entropyMatrix.length;
		int numberOfColumns = entropyMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns - 1); col++) {
				entropyMatrix[row][col] = calculateEntropyValue(row, col, grayscaleMatrix, grayscale9X9BlurMatrix);
			}
			for (int col = seamXValues[row] + 6; col < numberOfColumns - 1; col++) {
				entropyMatrix[row][col] = entropyMatrix[row][col + 1];
			}
		}
		removeLastColomn(entropyMatrix);
	}

	private void updateEdgeAndEntropyMatrix(int[] seamXValues) {
		int numberOfRows = edgeAndEntropyMatrix.length;
		int numberOfColumns = edgeAndEntropyMatrix[0].length;

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = Math.max(seamXValues[row] - 4, 0); col < Math.min(seamXValues[row] + 5,
					numberOfColumns - 1); col++) {
				edgeAndEntropyMatrix[row][col] = EDGES_WEIGHT * edgeMatrix[row][col]
						+ ENTROPY_WEIGHT * entropyMatrix[row][col];
			}
			for (int col = seamXValues[row] + 6; col < numberOfColumns - 1; col++) {
				edgeAndEntropyMatrix[row][col] = edgeAndEntropyMatrix[row][col + 1];
			}
		}
		removeLastColomn(edgeAndEntropyMatrix);
	}

	private static void removeLastColomn(double[][] array) {
		int numberOfRows = array.length;
		int numberOfColumns = array[0].length;
		for (int row = 0; row < numberOfRows; row++) {
			array[row] = Arrays.copyOfRange(array[row], 0, numberOfColumns - 1);
		}
	}

	private static void removeLastColomn(int[][][] array) {
		int numberOfRows = array.length;
		int numberOfColumns = array[0].length;
		for (int row = 0; row < numberOfRows; row++) {
			array[row] = Arrays.copyOfRange(array[row], 0, numberOfColumns - 1);
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
	
	
	public void enlargeImageHorizontallyByK(int k)
	{
		int[][] kMinSeams = SeamCarve.getKMinSeams(k,getEdgeAndEntropyMatrix(EnergyType.HoG));
		int[][][] resMatrix = new int[RGBMatrix.length][RGBMatrix[0].length+k][3];
		
		// insert new seams
		for(int i = 0; i < k; i++)
		{
			for(int row = 0; row < resMatrix.length; row++)
			{
				resMatrix[row][kMinSeams[i][row]] = new int[] {-1,-1,-1};
			}
		}
		
		// insert old matrix
		for(int row = 0; row < RGBMatrix.length; row++)
		{
			int offset = 0;
			for (int col = 0; col < RGBMatrix[0].length; col++)
			{
				while(resMatrix[row][col + offset][0] == -1)
				{
					offset++;
				} 
				resMatrix[row][col+offset] = RGBMatrix[row][col];
			}
		}
		
		// calculate averages for new seams:
		for(int row = 0; row < resMatrix.length; row++)
		{
			for (int col = 0; col<resMatrix[0].length; col++)
			{
				if(resMatrix[row][col][0]==-1)
				{
					if(col == 0)
					{
						resMatrix[row][col][0] = resMatrix[row][col+1][0]; 
						resMatrix[row][col][1] =  resMatrix[row][col+1][1];
						resMatrix[row][col][2] = resMatrix[row][col+1][2];
					} else if(col == resMatrix[0].length)
					{
						resMatrix[row][col][0] = resMatrix[row][col-1][0] ; 
						resMatrix[row][col][1] = resMatrix[row][col-1][1] ;
						resMatrix[row][col][2] = resMatrix[row][col-1][2] ;
					} else 
					{
						resMatrix[row][col][0] = (resMatrix[row][col-1][0] + resMatrix[row][col+1][0])/2; 
						resMatrix[row][col][1] = (resMatrix[row][col-1][1] + resMatrix[row][col+1][1])/2;
						resMatrix[row][col][2] = (resMatrix[row][col-1][2] + resMatrix[row][col+1][2])/2;
					}
				}
			}
		}
		RGBMatrix = resMatrix;
		
	}
	

}
