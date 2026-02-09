import javax.imageio.ImageIO;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
public class ImagePrep {

    // Resizes images to 64x64 with precision and not runtime preferred
    public Image sizer(BufferedImage image) {
        Image resizedImage = image.getScaledInstance(64, 64, Image.SCALE_SMOOTH);
        return resizedImage;
    }

    // Reads in images
    public Image readInImage(String name) {
        BufferedImage image = null;
        try {
            image = ImageIO.read(new File(name));
        } catch (IOException e) {}
        Image resizedImage = sizer(image);
        return resizedImage;
    }

    // Converts Image type to BufferedImage via Graphics
    public BufferedImage toBuffered(Image image) {
        BufferedImage bufferedImage = new BufferedImage(image.getWidth(null), image.getHeight(null), BufferedImage.TYPE_BYTE_GRAY);
        Graphics2D g = bufferedImage.createGraphics();
        g.drawImage(image, 0, 0, null);
        g.dispose();
        return bufferedImage;
    }

    //Converts the 64x64 image into an array 4096 grayscale values
    public double[] imageToData(String name) {
        Image image = readInImage(name);
        BufferedImage bimage = toBuffered(image);
        int[] pixels = bimage.getRGB(0, 0, bimage.getWidth(null), bimage.getHeight(null), null, 0, bimage.getWidth());
        double[] dpixels = new double[pixels.length];
        for (int i = 0; i < pixels.length; i++) {
            double gbyte = pixels[i] & 0xFF;
            dpixels[i] = gbyte / 255.0;
        }
        return dpixels;
    }
}
