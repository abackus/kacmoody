import org.jetbrains.annotations.Contract;

/**
 * Created by ABB on 5/6/2019.
 */
public class Math {
    // General math methods.
    @Contract(pure = true)
    static double bilinear(double[][] B, double[] x, double[] y, int n) {
        // Takes a n x n matrix B and two vectors x, y and computes the bilinear function x^T B y
        int o = 0;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                o += (x[i] * B[i][j] * y[j]);
            }
        }
        return o;
    }
}
