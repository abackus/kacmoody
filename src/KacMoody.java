import java.util.*;

public class KacMoody {
    // Class representing a Kac-Moody algebra.
    int dim; // The dimension of the maximal abelian subalgebra h.
    double[] functional; // The functional rho acting on h.
    double[][] Killing; // The Killing form B of the Kac-Moody algebra g.
    ArrayList<ArrayList<Root>> roots; // The positive roots of g, of each height.

    KacMoody(int dim, int height, double[] functional, double[][] Killing, ArrayList<Root> simpleRoots) {
        // TODO: Compute all these parameters from the Cartan matrix.
        this.dim = dim;
        this.functional = functional;
        this.Killing = Killing;
        this.roots = new ArrayList<>();
        this.roots.add(1, simpleRoots); // Simple roots have height 1.
        populateRoots(height);
    }

    private void populateRoots(int height) {
        //Generates positive roots up to given height.
        // TODO: This.
    }

    private double multiplicity(Root root) {
        //Computes the multiplicity of root using the Peterson recurrence formula.
        double p = peterson(root);
        // TODO: Remove the extra terms
        return p;
    }

    private double peterson(Root root) {
        // Given a root a computes the formally infinite sum c = sum_n m(a/n)/n
        if (root.mult != -1) {
            return root.mult;
        }

        double sum = 0;

        // TODO: Compute SUM

        double[] diff = new double[dim];
        for (int i=0;i<dim;i++) {
            diff[i] = 2 * functional[i] - root.coords[i];
        }

        return sum / Math.bilinear(Killing, root.coords, diff, dim);
    }
}
