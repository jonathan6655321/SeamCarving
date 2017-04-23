import java.awt.image.BufferedImage;

public class Matrix {
	public static void removeSeam(int[] seamXValues, int[][] matrix) {
		int numberOfRows = matrix.length;
		int numberOfColumns = matrix[0].length;
		int[][] newMatrix = new int[numberOfRows][numberOfColumns - 1];
		for (int row = 0; row < numberOfRows; row++) {
			System.arraycopy(matrix[row], 0, newMatrix[row], 0, seamXValues[row]);
			System.arraycopy(matrix[row], seamXValues[row] + 1, newMatrix[row], seamXValues[row],
					numberOfColumns - seamXValues[row] - 1);
			matrix[row] = newMatrix[row];
		}
	}

	public static void removeSeam(int[] seamXValues, double[][] matrix) {
		int numberOfRows = matrix.length;
		int numberOfColumns = matrix[0].length;
		double[][] newMatrix = new double[numberOfRows][numberOfColumns - 1];
		for (int row = 0; row < numberOfRows; row++) {
			System.arraycopy(matrix[row], 0, newMatrix[row], 0, seamXValues[row]);
			System.arraycopy(matrix[row], seamXValues[row] + 1, newMatrix[row], seamXValues[row],
					numberOfColumns - seamXValues[row] - 1);
			matrix[row] = newMatrix[row];
		}
	}

	public static double[][] rotateArray(double[][] array) {
		int numberOfRows = array.length;
		int numberOfColumns = array[0].length;
		double[][] rotatedArray = new double[numberOfColumns][numberOfRows];
		for (int row = 0; row < numberOfRows; row++) {
			for (int col = 0; col < numberOfColumns; col++) {
				rotatedArray[col][numberOfRows - row - 1] = array[row][col];
			}
		}
		return rotatedArray;
	}

	public static int[][] rotateArray(int[][] array) {
		int numberOfRows = array.length;
		int numberOfColumns = array[0].length;
		int[][] rotatedArray = new int[numberOfColumns][numberOfRows];
		for (int row = 0; row < numberOfRows; row++) {
			for (int col = 0; col < numberOfColumns; col++) {
				rotatedArray[col][numberOfRows - row - 1] = array[row][col];
			}
		}
		return rotatedArray;
	}

	public static BufferedImage createBufferImageFromRGBMatrix(int[][] RGBMatrix) {
		int numberOfRows = RGBMatrix.length;
		int numberOfColumns = RGBMatrix[0].length;
		BufferedImage image = new BufferedImage(numberOfColumns, numberOfRows, BufferedImage.TYPE_INT_RGB);

		for (int row = 0; row < numberOfRows; row++) {
			for (int col = 0; col < numberOfColumns; col++) {
				image.setRGB(col, row, RGBMatrix[row][col]);
			}
		}
		return image;
	}

	public static BufferedImage createBufferImageFromIntMatrix(int[][] matrix) {
		int height = matrix.length;
		int width = matrix[0].length;
		int[][] newRGBMatrix = new int[height][width];
		double maxValue = 0;
		for (int[] dArr : matrix) {
			for (int d : dArr) {
				maxValue = Math.max(maxValue, d);
			}
		}
		System.out.println(maxValue);
		maxValue = ((double) 250) / maxValue;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				newRGBMatrix[i][j] = (int) (maxValue * matrix[i][j]);
				newRGBMatrix[i][j] += (newRGBMatrix[i][j] << 16) + (newRGBMatrix[i][j] << 8);
			}
		}
		return Matrix.createBufferImageFromRGBMatrix(newRGBMatrix);
	}

	public static BufferedImage createBufferImageFromDoubleMatrix(double[][] matrix) {
		int height = matrix.length;
		int width = matrix[0].length;
		int[][] newRGBMatrix = new int[height][width];
		double maxValue = 0;
		for (double[] dArr : matrix) {
			for (double d : dArr) {
				maxValue = Math.max(maxValue, d);
			}
		}
		System.out.println(maxValue);
		maxValue = ((double) 250) / maxValue;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				newRGBMatrix[i][j] = (int) (maxValue * matrix[i][j]);
				newRGBMatrix[i][j] += (newRGBMatrix[i][j] << 16) + (newRGBMatrix[i][j] << 8);
			}
		}
		return Matrix.createBufferImageFromRGBMatrix(newRGBMatrix);
	}

	public static double[][] clone(double[][] original) {
		if (original.length == 0)
			return new double[0][0];
		double[][] newMat = new double[original.length][original[0].length];
		for (int i = 0; i < original.length; i++)
			System.arraycopy(original[i], 0, newMat[i], 0, original[i].length);
		return newMat;
	}

	public static int[][] clone(int[][] original) {
		if (original.length == 0)
			return new int[0][0];
		int[][] newMat = new int[original.length][original[0].length];
		for (int i = 0; i < original.length; i++)
			System.arraycopy(original[i], 0, newMat[i], 0, original[i].length);
		return newMat;
	}
}
