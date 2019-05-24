# Example of use:
# CartanMatrix([[2,-3],[-3,2]]).rm(height=25)


def sub_parts(c):
    """
    For a list c = [1,2,0], return
                  [[0,0,0],
                   [1,0,0],
                   [0,1,0],
                   [1,1,0],
                   [0,2,0]]
    """
    #TODO: It might be possible to make this faster with memoization -- but
    # not much faster, since the set of subparts is *very large*.
    # Better to just use the Weyl group tbh.
    if len(c) == 1:
        return [[x] for x in range(c[0])]

    if c[0] == 0:
        return [[0] + y for y in sub_parts(c[1:])]

    d = sub_parts(c[1:])
    g = sub_parts([c[0]-1]+c[1:])+[[c[0]]+y for y in sub_parts(c[1:])]+[[c[0]-1]+c[1:]]
    return g

def common_divisors(l):
    # Returns the common divisors of a list l of integers >= 0.
    return divisors(gcd(l))

def all_lists(l,s):
    #return all lists with length l and sum s
    if l == 1:
        return [[s]]
    to_return = []
    for i in range(s+1):
        to_return = to_return + [[i]+x for x in all_lists(l-1,s-i)]
    return to_return

def rm(self, height=10):
    # Returns the root multiplicities of a Cartan matrix self
    # Recall that if beta is a root, then mult(beta) = mult(-beta)
    # So we only have to check positive roots.

    # Peterson's recurrence formula:
    # For beta a positive root,
    # B(beta, 2rho - beta)c_beta = sum_{alpha+gamma=beta} B(alpha, gamma)c_alphac_gamma
    # where B is the Killing form, rho is the is the functional that sends the coroot
    # hat alpha_i to a_ii/2, and c(beta) = sum_n mult(beta/n)/n.

    self.multiplicities = {}
    self.c = {}
    self.c[tuple([0, 0, 0])] = 0
    mat = self.symmetrized_matrix() # TODO: This is probably a BUG if mat is not symmetric. Fix it.
    dim = mat.nrows()

    def two_rho(beta):
        # The functional rho which appears in Peterson's recurrence formula
        # TODO: This is probably a BUG if mat is not symmetric. Fix it.
        return -2*sum(beta)

    def B(a, b):
        # The Killing form of mat
        # TODO: This is probably a BUG if mat is not symmetric. Fix it.
        return a * mat * b

    # Initialize the dict of root multiplicities so we can start the induction.
    for root in all_lists(dim, 1):
        # Iterate over all positive simple roots of mat.
        self.multiplicities[tuple(root)] = 1
        self.c[tuple(root)] = 1

    # Now do the inductive step.
    # TODO: Get rid of all these tuple casts
    for h in range(2, height + 1):
        for root in all_lists(dim, h):
            # Iterate over all positive elements of the root lattice of mat,
            # with height ( = ell^1 norm) h >= 2.
            RHS = 0 # Right-hand side of Peterson's recurrence formula
            root = vector(root)
            for subroot in sub_parts(list(root)):
                # Sum over subroots
                # TODO: We don't actually need to sum over the subroots because the
                # action of the Weyl group preserves multiplicities, and hence c.
                # So this could be greatly sped up.
                subroot = vector(subroot)
                if subroot == 0:
                    continue
                #print('\t' + str(subroot))

                RHS = RHS + B(subroot, root - subroot)*self.c[tuple(subroot)]*self.c[tuple(root - subroot)]

            # The factor which appears in the left-hand side of Peterson's recurrence formula
            # TODO: This is probably a BUG if mat is not symmetric. Fix it.
            LHSFactor = B(root, root) + two_rho(root)

            # Compute c(root) - m(root)
            cRoot = 0
            for n in divisors(gcd(root)):
                if n != 1:
                    cRoot = cRoot + self.multiplicities[tuple(root/n)]/n

            if RHS == 0 and LHSFactor == 0:
                # Not a root
                self.multiplicities[tuple(root)] = 0
                self.c[tuple(root)] = cRoot
                continue

            self.multiplicities[tuple(root)] = RHS/LHSFactor - cRoot
            self.c[tuple(root)] = RHS/LHSFactor

    # Clean the root lattice, removing zero entries and replacing floats with ints
    notRoots = [] # Can't change the size of a dict during iteration
    for r in self.multiplicities:
        if self.multiplicities[r] == 0:
            notRoots.append(r)
    for r in notRoots:
        self.multiplicities.pop(r)

    return self.multiplicities


CartanMatrix.rm = rm



def CartanMatrixGenerator(dim):
    # Iterates over Cartan matrices of dimension dim
    # Create 2I_dim
    matrix = [[0 for i in range(dim)] for i in range(dim)]
    for i in range(dim):
        matrix[i][i] = 2
    order = 0 # By "order" of I mean the sum of -a_{ij}, i \neq j
    while true:
        # It might not be necessary to check *every* matrix, by checking
        # if a matrix is a direct sum. TODO this
        order = order + 1
        for vector in all_lists(dim * (dim - 1), order):
            k = 0
            for i in range(dim):
                for j in range(dim):
                    if i != j:
                        matrix[i][j] = -vector[k]
                        k = k + 1
            symmetricZeroes = true
            for i in range(dim):
                for j in range(dim): # This could be optimized...
                    if (matrix[i][j] == 0) and (matrix[j][i] != 0):
                        symmetricZeroes = false
            if symmetricZeroes:
                yield CartanMatrix(matrix)



def numerology(dim, count, h=10):
    # Computes root multiplicities for dimension dim with height h, until count Kac-Moody algebras have been studied.
    with open('data.csv', 'a') as file:
        gen = CartanMatrixGenerator(dim)
        for i in range(count):
            kacmoody = next(gen)
            if kacmoody.is_symmetrizable():
                roots = kacmoody.rm(height=h)
                file.write(str(kacmoody) + "\n")
                print(str(kacmoody) + "\n")
                for root in roots:
                    try:
                        file.write(str(root) + ", " + str(roots[root]) + "\n")
                    except ZeroDivisionError as (errno, errstr):
                        print "ZeroDivisionError({0}): {1}".format(errno, strerror)
                file.write("\n\n")
