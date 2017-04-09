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

import javax.imageio.ImageIO;
import javax.swing.JComponent;
import javax.swing.JFrame;

import java.awt.Component;
import java.awt.image.BufferedImage;
import java.io.IOException;
import javax.imageio.ImageIO;

public class SeamCarve {
	private SeamImage originalImage;
	public static void main(String[] args) {
		SeamCarve SC = new SeamCarve();
		SC.resizePicture("C:\\Users\\kessi\\Documents\\Tau\\Graphic\\SeamCurve\\bin\\reef-view-hotel-5.jpg", 0, 0,
				EnergyType.HoG, "lol");
	}

	private void resizePicture(String imageFilename, int numCol, int numRow, EnergyType eType,
			String outputImageFilename) {
		originalImage = new SeamImage(imageFilename);
		displayImage(originalImage.getOriginalImage());		
		
	}


	public void displayImage(BufferedImage image) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				ImageFrame frame = new ImageFrame(image);
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.setVisible(true);

			}
		});
	}

	

	class ImageFrame extends JFrame {

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
	
	public int[][] calculateMinSeamsMatrix(int[][] edgeMatrix)
	{
		int[][] minSeamsMatrix = matrixCopy(edgeMatrix);
		int edgeMatrixWidth = edgeMatrix[0].length;
		int edgeMatrixHeight = edgeMatrix.length;
		
		for (int i=1; i<edgeMatrixHeight; i++)
		{
			for (int j=0; j<edgeMatrixWidth; j++)
			{
				LinkedList<Integer> minPathBeforIJ = new LinkedList<>();
				minPathBeforIJ.add(minSeamsMatrix[i-1][j]);
				if(j!=0){
					minPathBeforIJ.add(minSeamsMatrix[i-1][j-1]);
				}
				if(j!=edgeMatrixWidth-1)
				{
					minPathBeforIJ.add(minSeamsMatrix[i-1][j+1]);
				}
				minSeamsMatrix[i][j] = Collections.min(minPathBeforIJ) + edgeMatrix[i][j];
			}
		}
		return minSeamsMatrix;
	}
	
	public int[] getMinSeam(int[][] minSeamsMatrix)
	{
		int numCol = minSeamsMatrix[0].length;
		int numRows = minSeamsMatrix.length;
		int[] seam = new int[numRows]; 
		seam[numRows-1] = minElementsIndex(minSeamsMatrix[numRows-1]);
		
		for(int i=numRows-2; i>=0; i--)
		{
			int prevIndex = seam[i+1];
			
			int minVal = minSeamsMatrix[i][prevIndex]; // directly above
			int minIndex = prevIndex;
			
			if(prevIndex!=0 && minVal > minSeamsMatrix[i][prevIndex-1])
			{
				minVal = minSeamsMatrix[i][prevIndex-1];
				minIndex--;
			}
			if(prevIndex!=numCol-1 && minVal > minSeamsMatrix[i][prevIndex+1])
			{
				minVal = minSeamsMatrix[i][prevIndex+1];
				minIndex = prevIndex + 1;
			}
			
			seam[i] = minIndex;
		}
		return seam;
	}
	
	public int minElementsIndex(int[] arr)
	{
		int min = arr[0];
		int j = 0;
		
		for (int i=0; i<arr.length; i++)
		{
			if(arr[i] < min)
			{
				min = arr[i];
				j = i;
			}		
		}
		
		return j;
	}
	
	
	public int[][] matrixCopy(int[][] original)
	{
		int [][] newMat = new int[original.length][];
		for(int i = 0; i < original.length; i++)
			newMat[i] = original[i].clone();
		return newMat;
	}
}
