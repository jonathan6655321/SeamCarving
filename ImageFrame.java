import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Insets;
import java.awt.image.BufferedImage;

import javax.swing.JComponent;
import javax.swing.JFrame;

class ImageFrame extends JFrame {
	private static final long serialVersionUID = 1L;
	public ImageComponent IC;

	public static void displayImage(SeamImage image) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				ImageFrame frame = new ImageFrame(image.getImage());
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.setVisible(true);

			}
		});
	}

	private ImageFrame(BufferedImage image1) {
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

			g.drawImage(image, 0, 0, this);
		}

	}
}