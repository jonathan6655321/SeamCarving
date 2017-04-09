import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.io.File;
import java.io.IOException;
import java.util.Observable;
import java.util.Observer;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
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
}
