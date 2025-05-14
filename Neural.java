import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class Neural {
    long seed = System.nanoTime();
    Random random = new Random(seed);
    int[] n = {4096, 128, 64, 1};
    double[][] W1 = randn(n[1], n[0]);
    double[][] W2 = randn(n[2], n[1]);
    double[][] W3 = randn(n[3], n[2]);
    double[][] b1 = randn(n[1], 1);
    double[][] b2 = randn(n[2], 1);
    double[][] b3 = randn(n[3], 1);
    double[][] X = new double[4096][200];
    double[][] A0 = new double[4096][200];
    double[][] y = new double[200][1];
    double[][] Y = transpose(y);
    double[] mean = new double [4096];
    double[] stdev = new double[4096];
    class matrixTrio {
        double[][] W;
        double[][] b;
        double[][] A;
        matrixTrio(double[][] W, double[][] b, double[][] A)  {
            this.W = W;
            this.b = b;
            this.A = A;
        }
    }
    class matrixDuo {
        double[][] W;
        double[][] b;
        matrixDuo(double[][] W, double[][] b) {
            this.W = W;
            this.b = b;
        }
    }
    class AMatricies {
        double[][] A3;
        double[][] A2;
        double[][] A1;
        double[][] A0;
        AMatricies(double[][] A3, double[][] A2, double[][] A1, double[][] A0) {
            this.A3 = A3;
            this.A2 = A2;
            this.A1 = A1;
            this.A0 = A0;
        }
    }
    public void standardA0() {
        this.A0 = standardize(this.X);
    }
    public double[][] setY(double[][] m) {
        for (int i = 0; i < 100; i++) {
            m[i][0] = 1;
        }
        for (int i = 100; i < 200; i++) {
            m[i][0] = 0;
        }
        return m;
    }
    public double[][] addToX(double[] insert, int j) {
            for (int i = 0; i < insert.length; i++) {
                X[i][j] = insert[i];
            }
        return X;
    }
    public double[][] randn(int rows, int columns) {
        double[][] m = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                m[i][j] = random.nextGaussian();
            }
        }
        return m;
    }
    public double[][] transpose(double[][] m) {
        double[][] mtrans = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                mtrans[j][i] = m[i][j];
            }
        }
        return mtrans;
    }
    public int[][] transpose(int[][] m) {
        int[][] mtrans = new int[m[0].length][m.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                mtrans[j][i] = m[i][j];
            }
        }
        return mtrans;
    }
    public String shape(double[][] m) {
        return "(" + m.length + ", " + m[0].length + ")";
    }
    public String shape(int[][] m) {
        return "(" + m.length + ", " + m[0].length + ")";
    }
    public double[][] sigmoid(double[][] m) {
        double[][] result = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                result[i][j] = 1 / (1 + Math.exp(-1 * m[i][j]));
            }
        }
        return result;
    }
    public double[][] ln(double[][] m) {
        double[][] result = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                result[i][j] = Math.log(m[i][j]);
            }
        }
        return result;
    }
    public double[][] oneMinus(double[][] m) {
        double[][] result = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                result[i][j] = 1 - m[i][j];
            }
        }
        return result;
    }
    public double[][] negative(double[][] m) {
        double[][] result = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                result[i][j] = -1 * m[i][j];
            }
        }
        return result;
    }
    public double[][] multiply(double[][] m1, double[][] m2) {
        double[][] result = new double[m1.length][m2[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m2[0].length; j++) {
                double res = 0;
                for (int k = 0; k < m1[0].length; k++) {
                    res += m1[i][k] * m2[k][j];
                }
                result[i][j] = res;
            }
        }
        return result;
    }
    public double[][] add(double[][] m1, double [][] m2) {
        int columns = Math.max(m1[0].length, m2[0].length);
        double[][] result = new double[m1.length][columns];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < columns; j++) {
            result[i][j] = m1[i][j] + m2[i][0];
            }
        }
        return result;
    }
    public double[][] sum(double[][] m) {
        double[][] result = new double[m.length][1];
        for (int i = 0; i < m.length; i++) {
            double sum = 0;
            for (int j = 0; j < m[0].length; j++) {
                sum += m[i][j];
            }
            result[i][0] = sum;
        }
        return result;
    }
    public double cost(double[][] y_hat, double[][] y) {
        int m = y_hat[0].length;
        double losses = 0;
        for (int i = 0; i < m; i++) {
            double yi = y[0][i];
            double y_hati = y_hat[0][i];
            losses += - ((yi * Math.log(y_hati)) + (1 - yi) * Math.log(1 - y_hati));
        }
        double summed_losses = losses / m;
        return summed_losses;
    }
    public double[][] standardize(double[][] m) {
        double[][] result = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            double sum = 0;
            double stdevsum = 0;
            int j = 0;
            for (j = 0; j < m[i].length; j++) {
                sum += m[i][j];
            }
            mean[i] = sum / j;
            for (int k = 0; k < m[i].length; k++) {
                stdevsum += Math.pow(m[i][k] - mean[i], 2);
            }
            stdev[i] = Math.sqrt(stdevsum / j);
            stdev[i] = stdev[i] == 0 ? 1e-8 : stdev[i];
            for (int l = 0; l < m[i].length; l++) {
                result[i][l] = (m[i][l] - mean[i]) / stdev[i];
            } 
        }
        return result;
    }
    public matrixTrio backpropLayer3(double[][] y_hat, double[][] Y, int m, double[][] A2, double[][] W3) {
        double[][] A3 = y_hat;
        m = A3[0].length;
        double[][] dC_dZ3 = new double[1][m];
        for (int i = 0; i < m; i++) {
            dC_dZ3[0][i] = (A3[0][i] - Y[0][i]) / m;
        }
        double[][] dZ3_dW3 = A2;
        double[][] dC_dW3 = multiply(dC_dZ3, transpose(dZ3_dW3));
        double[][] dC_db3 = sum(dC_dZ3); 
        double[][] dZ3_dA2 = W3;
        double[][] dC_dA2 = multiply(transpose(W3), dC_dZ3);
        matrixTrio trio = new matrixTrio(dC_dW3, dC_db3, dC_dA2);
        return trio;
    }
    public matrixTrio backpropLayer2(double[][]propagator_dC_dA2, double[][]A1, double[][]A2, double[][]W2) {
        double[][] dA2_dZ2 = new double[A2.length][A2[0].length];
        for (int i = 0; i < A2.length; i++) {
            for (int j = 0; j < A2[0].length; j++) {
                dA2_dZ2[i][j] = A2[i][j] * (1 - A2[i][j]);
            }
        }
        double[][] dC_dZ2  = new double[A2.length][A2[0].length];
        for (int i = 0; i < A2.length; i++) {
            for (int j = 0; j < A2[0].length; j++) {
            dC_dZ2[i][j] = propagator_dC_dA2[i][j] * dA2_dZ2[i][j];
            }
        }
        double[][] dZ2_dW2 = A1;
        double[][] dC_dW2 = multiply(dC_dZ2, transpose(dZ2_dW2));
        double[][] dC_db2 = sum(dC_dZ2);
        double[][] dZ2_dA1 = W2;
        double[][] dC_dA1 = multiply(transpose(W2), dC_dZ2);
        matrixTrio trio = new matrixTrio(dC_dW2, dC_db2, dC_dA1);
        return trio;
    }
    public matrixDuo backpropLayer1(double[][] propagator_dC_dA1, double[][] A1, double[][] A0, double[][] W1) {
        double[][] dA1_dZ1 = new double[A1.length][A1[0].length];
        for (int i = 0; i < A1.length; i++) {
            for (int j = 0; j < A1[0].length; j++) {
                dA1_dZ1[i][j] = A1[i][j] * (1 - A1[i][j]);
            }
        }
        double[][] dC_dZ1 = new double[A1.length][A1[0].length];
        for (int i = 0; i < A1.length; i++) {
            for (int j = 0; j < A1[0].length; j++) {
                dC_dZ1[i][j] = propagator_dC_dA1[i][j] * dA1_dZ1[i][j];
            }
        }
        double[][] dZ1_dW1 = A0;
        double[][] dC_dW1 = multiply(dC_dZ1, transpose(dZ1_dW1));
        double[][] dC_db1 = sum(dC_dZ1);
        matrixDuo duo = new matrixDuo(dC_dW1, dC_db1);
        return duo;
    }
    public AMatricies feedForward(double[][] A0) {
        double[][] Z1 = add(multiply(W1, A0), b1);
        double[][] A1 = sigmoid(Z1);
        double[][] Z2 = add(multiply(W2, A1), b2);
        double[][] A2 = sigmoid(Z2);
        double[][] Z3 = add(multiply(W3, A2), b3);
        double[][] A3 = sigmoid(Z3); 
        AMatricies As = new AMatricies(A3, A2, A1, A0);
        return As;
    }
    public double[][] gradientUpdates(double[][] m, double alpha, double[][] n) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                m[i][j] = m[i][j] - (alpha * n[i][j]);
            }
        }
        return m;
    }
    public ArrayList<Double> train() {
        int epochs = 100;
        double alpha = 1;
        ArrayList<Double> costs = new ArrayList<>();
        for (int i = 0; i < epochs; i++) {
            int m = A0[0].length;
            AMatricies As = feedForward(A0);
            double[][] y_hat = As.A3;
            double error = cost(y_hat, Y);
            costs.add(error);
            matrixTrio trio1 = backpropLayer3(y_hat, Y, m, As.A2, W3);
            double[][] dC_dW3 = trio1.W;
            double[][] dC_db3 = trio1.b;
            double[][] dC_A2 = trio1.A;
            matrixTrio trio2 = backpropLayer2(dC_A2, As.A1, As.A2, W2);
            double[][] dC_dW2 = trio2.W;
            double[][] dC_db2 = trio2.b;
            double[][] dC_A1 = trio2.A;
            matrixDuo duo1 = backpropLayer1(dC_A1, As.A1, As.A0, W1);
            double[][] dC_dW1 = duo1.W;
            double[][] dC_db1 = duo1.b;
            W3 = gradientUpdates(W3, alpha, dC_dW3);
            W2 = gradientUpdates(W2, alpha, dC_dW2);
            W1 = gradientUpdates(W1, alpha, dC_dW1);
            b3 = gradientUpdates(b3, alpha, dC_db3);
            b2 = gradientUpdates(b2, alpha, dC_db2);
            b1 = gradientUpdates(b1, alpha, dC_db1);

                System.out.printf("epoch %3d  cost %.5f%n", i, error);
        }
        return costs;
    }
    public static void main(String[] args) {
        Neural net = new Neural();
        ImagePrep imagePrep = new ImagePrep();
        for (int i = 1; i <= 100; i++) {
            String name = "C:\\Users\\jerom\\Downloads\\burger" + i + ".jpg";
            net.addToX(imagePrep.imageToData(name), i - 1);
        }
        for (int i = 1; i <=100; i++) {
            String name = "C:\\Users\\jerom\\Downloads\\nonburger" + i + ".jpg";
            net.addToX(imagePrep.imageToData(name), 100 + (i - 1));
        }
        net.standardA0();
        net.y = net.setY(net.y);
        net.Y = net.transpose(net.y);
        net.train();

        // forward pass with the trained weights
        AMatricies As = net.feedForward(net.A0);
        double[][] A3 = As.A3;            // shape (1, m = 10)
        System.out.println(Arrays.toString(A3[0]));

        for (int i = 0; i < 100; i++) {
            String name = "C:\\Users\\jerom\\Downloads\\testburger" + (i + 1) + ".jpg";
            double[][] testm = new double[4096][1];
            double[] imagea = imagePrep.imageToData(name);
            for (int j = 0; j < testm.length; j++) {
                testm[j][0] = (imagea[j] - net.mean[j]) / net.stdev[j];
            }
            AMatricies testAs = net.feedForward(testm);
            double[][] A3test = testAs.A3;
            if (A3test[0][0] < 0.5) {
                System.out.println("Image " + (i + 1) + " is not a burger. A3test[0][0] is " + A3test[0][0]);
            }
            else {
                System.out.println("Image " + (i + 1) +" is a burger. A3test[0][0] is " + A3test[0][0]);
            }
        }
    }
}
