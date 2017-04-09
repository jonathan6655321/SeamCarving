import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Observable;
import java.util.Observer;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
import java.awt.List;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;

import javax.imageio.ImageIO;
import javax.swing.JComponent;
import javax.swing.JFrame;

import java.awt.Component;
import java.awt.image.BufferedImage;
import java.io.IOException;
import javax.imageio.ImageIO;
/* TODO:
 * enlarage bigger then original*2;
 */
public class SeamCarve {
	private SeamImage originalImage;

	public static void main(String[] args) {
		SeamCarve SC = new SeamCarve();
		SC.resizePicture(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]), EnergyType.HoG, args[4]);
	}

	private void resizePicture(String imageFilename, int numCol, int numRow, EnergyType eType,
			String outputImageFilename) {
		originalImage = new SeamImage(imageFilename);

		if (numCol > originalImage.getWidth()) {
			originalImage.enlargeImageHorizontallyByK(getKMinSeams(numCol - originalImage.getWidth(), originalImage));
		}
		if (numRow > originalImage.getHeight()) {
			originalImage.rotate90right();
			originalImage.enlargeImageHorizontallyByK(getKMinSeams(numRow - originalImage.getWidth(), originalImage));
			originalImage.rotate90right();
			originalImage.rotate90right();
			originalImage.rotate90right();
		}
		if (numCol < originalImage.getWidth()) {
			int numberOfSeamToRemove =originalImage.getWidth() - numCol;
			for (int i = 0; i <numberOfSeamToRemove ; i++) {
				originalImage.removeVerticalSeam(
						getMinSeam(calculateMinSeamsMatrix(originalImage.getEdgeAndEntropyMatrix(EnergyType.HoG))));
			}
		}

		if (numRow < originalImage.getHeight()) {
			originalImage.rotate90right();
			int numberOfSeamToRemove =originalImage.getWidth() - numRow;
			for (int i = 0; i < numberOfSeamToRemove; i++) {
				originalImage.removeVerticalSeam(
						getMinSeam(calculateMinSeamsMatrix(originalImage.getEdgeAndEntropyMatrix(EnergyType.HoG))));
			}
			originalImage.rotate90right();
			originalImage.rotate90right();
			originalImage.rotate90right();
		}

		displayImage(originalImage);

	}

	public void displayImage(SeamImage image) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				ImageFrame frame = new ImageFrame(image.getOriginalImage());
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.setVisible(true);
				// ImageComponent IC = frame.getImageComponent();
				//
				// // frame.addComponentListener(new ComponentAdapter()
				// // {
				// // public void componentResized(ComponentEvent evt) {
				// // Component c = (Component)evt.getSource();
				// // //........
				// // }
				// // });
				// // for (int i = 0; i <300; i++) {
				// //// try {
				// //// Thread.sleep(250);
				// //// } catch (InterruptedException e) {
				// //// // TODO Auto-generated catch block
				// //// e.printStackTrace();
				// //// }
				// // int[] seam =
				// //
				// getMinSeam(calculateMinSeamsMatrix(image.getEdgeAndEntropyMatrix(EnergyType.HoG)));
				// // image.removeVerticalSeam(seam);
				// // }
				// // IC.setImage(image.getOriginalImage());

			}
		});
	}

	class ImageFrame extends JFrame {
		public ImageComponent IC;

		public ImageFrame(BufferedImage image1) {
			JFrame temp = new JFrame();
			temp.pack();
			Insets insets = temp.getInsets();
			temp = null;
			setSize(new Dimension(insets.left + insets.right + image1.getWidth(),
					insets.top + insets.bottom + image1.getHeight()));
			this.setVisible(true);
			this.setResizable(false);

			ImageComponent component = new ImageComponent(image1);
			add(component);
			IC = component;

		}

		public ImageComponent getImageComponent() {
			return IC;
		}
	}

	class ImageComponent extends JComponent {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private Image image;

		public ImageComponent(Image image1) {
			image = image1;
		}

		public void setImage(Image image1) {
			image = image1;
		}

		public void paintComponent(Graphics g) {
			if (image == null)
				return;
			int imageWidth = image.getWidth(this);
			int imageHeight = image.getHeight(this);

			g.drawImage(image, 0, 0, this);

			// for (int i = 0; i * imageWidth <= getWidth(); i++)
			// for (int j = 0; j * imageHeight <= getHeight(); j++)
			// if (i + j > 0)
			// g.copyArea(0, 0, imageWidth, imageHeight, i * imageWidth, j *
			// imageHeight);
		}

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

	public static int[] getMinSeam(double[][] minSeamsMatrix) {
		int numCol = minSeamsMatrix[0].length;
		int numRows = minSeamsMatrix.length;
		int[] seam = new int[numRows];
		seam[numRows - 1] = minElementsIndex(minSeamsMatrix[numRows - 1]);

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

	public static int[][] getKMinSeams(int k, SeamImage originalImage) {
		int numRows = originalImage.getEdgeAndEntropyMatrix(EnergyType.HoG).length;
		int numCols = originalImage.getEdgeAndEntropyMatrix(EnergyType.HoG)[0].length;
		int[][] kMinSeams = new int[k][numRows];
		SeamImage image = new SeamImage(originalImage.getOriginalImage());

		boolean[][] isInSeam = new boolean[numRows][numCols];

		for (int i = 0; i < k; i++) {
			kMinSeams[i] = getMinSeam(calculateMinSeamsMatrix(image.getEdgeAndEntropyMatrix(EnergyType.HoG)));
			image.removeVerticalSeam(kMinSeams[i]);
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

		return calcTrueInsertedIndex(kMinSeams);
	}

	public static int[][] calcTrueInsertedIndex(int[][] kMinSeams) {
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

	public static int minElementsIndex(double[] arr) {
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
