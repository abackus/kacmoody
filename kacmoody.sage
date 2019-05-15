# Example of use:
# CartanMatrix([[2,-3],[-3,2]]).rm(depth=25)


def sub_parts(c):
    """
    For a list c = [1,2,0], return
                  [[0,0,0],
                   [1,0,0],
                   [0,1,0],
                   [1,1,0],
                   [0,2,0]]
    """
    if len(c) == 1:
        return [[x] for x in range(c[0])]

    if c[0] == 0:
        return [[0] + y for y in sub_parts(c[1:])]

    d = sub_parts(c[1:])
    g = sub_parts([c[0]-1]+c[1:])+[[c[0]]+y for y in sub_parts(c[1:])]+[[c[0]-1]+c[1:]]
    return g

def common_divisors(l):
    # Returns the common divisors of a list l of integers >= 0.
    l = [a for a in l if a != 0]
    d = [set(e.divisors()) for e in l]
    f = d[0]
    for i in d[1:]:
        f = f.intersection(i)
    return f

def all_lists(l,s):
    #return all lists with length l and sum s
    if l == 1:
        return [[s]]
    to_return = []
    for i in range(s+1):
        to_return = to_return + [[i]+x for x in all_lists(l-1,s-i)]
    return to_return



def rm(self,depth=5):
    # Returns the root multiplicities of a Cartan matrix self
    # Recall that if beta is a root, then mult(beta) = mult(-beta)
    # So we only have to check positive roots.
    R = self.root_system()
    L = R.root_lattice()
    S = L.simple_roots()
    self.multiplicities = {}
    D = self.symmetrized_matrix().dense_matrix()

    # Memoization of c
    cMemo = {}

    # Peterson's recurrence formula:
    # For beta a positive root,
    # B(beta, 2rho - beta)c_beta = sum_{alpha+gamma=beta} B(alpha, gamma)c_alphac_gamma
    # where B is the Killing form, rho is the is the functional that sends the coroot
    # hat alpha_i to a_ii/2, and c_beta = sum_n mult(beta/n)/n.

    def c(beta):
        if beta in cMemo:
            return cMemo[beta]
        # The sum in Peterson's recurrence formula is formally infinite, but
        # since only when div|beta is the value nonzero, we only have to sum then.
        divisors = common_divisors(beta.dense_coefficient_list())
        c = 0
        for div in divisors:
            if beta/div in self.multiplicities:
                c = c + (float(1)/float(div))*self.multiplicities[beta/div]
        cMemo.update({beta: c})
        return c

    def two_rho(beta):
        return -2*sum(beta.dense_coefficient_list())

    def c_to_root(coeff):
        # Given coefficients, return the root they describe.
        #IT IS 0 INDEXED FOR INFINITE, 1 INDEXED FOR FINITE
        #THIS IS A BUG IN SAGE
        if 0 in S.keys():
            s = sum([S[int(x)]*int(coeff[x]) for x in range(len(coeff))])
        else:
            s = sum([S[int(x+1)]*int(coeff[x]) for x in range(len(coeff))])
        return s

    # Generates the root lattice up to the given depth.
    #Actually, it is SO much easier to induct on height...
    all_combs = []
    for l in [all_lists(len(S),x) for x in range(1,depth+1)]:
        all_combs = all_combs + l
    roots = [c_to_root(comb) for comb in all_combs]

    for r in roots:
        rd = r.dense_coefficient_list()
        if sum(rd) == 1: # If r is a trivial root
            self.multiplicities[r] = 1
        else:
            a_o = 0
            for bp in [c_to_root(l) for l in sub_parts(rd) if sum(l) != 0]:
                # Sum over subroots. Could this be memoized?
                bpp = r-bp
                bpd = bp.dense_coefficient_list()
                bppd = bpp.dense_coefficient_list()
                a_o = a_o + vector(bppd)*D*vector(bpd)*c(bp)*c(bpp)
            b = vector((r).dense_coefficient_list())*D*vector((r).dense_coefficient_list()) + two_rho(r)

            c_o = 0
            for n in common_divisors(rd):
                if n != 1:
                    c_o = c_o + float(1)/float(n)*self.multiplicities[c_to_root(list(vector(rd)/n))]


            if a_o == 0 and b == 0:
                self.multiplicities[r] = 0
                continue
            self.multiplicities[r] = a_o/b-c_o

    # Clean the root lattice, removing zero entries and replacing floats with ints
    notRoots = [] # Can't change the size of a dict during iteration
    for r in self.multiplicities:
        self.multiplicities[r] = int(self.multiplicities[r])
        if self.multiplicities[r] == 0:
            notRoots.append(r)
    for r in notRoots:
        self.multiplicities.pop(r)
    return self.multiplicities

CartanMatrix.rm = rm

def CartanMatrixGenerator(dim):
    # Iterates over Cartan matrices of dimension dim
    # TODO check: this is surjective
    matrix = []
    for i in range(dim):
        row = []
        for j in range(dim):
            if i == j:
                row.append(2)
            else:
                row.append(0)
        matrix.append(row) # Construct the trivial Cartan matrix 2I
    yield CartanMatrix(matrix)

    #TODO make this work in the nontrivial cases
    # But we're just doing this for 2x2 symmetric matrices for now
    while true:
        matrix[0][1] = matrix[0][1] - 1
        matrix[1][0] = matrix[1][0] - 1
        yield CartanMatrix(matrix)


def numerology(dim, count, dep):
    # Computes root multiplicities for dimension dim with depth dep, until count Kac-Moody algebras have been studied.
    with open('data.csv', 'a') as file:
        gen = CartanMatrixGenerator(dim)
        for i in range(count):
            kacmoody = next(gen)
            if kacmoody.is_symmetrizable():
                roots = kacmoody.rm(depth=dep)
            file.write(str(kacmoody) + "\n") #TODO debug
            for root in roots:
                try:
                    file.write(str(root) + ", " + str(roots[root]) + "\n")
                except ZeroDivisionError as (errno, errstr):
                    print "ZeroDivisionError({0}): {1}".format(errno, strerror)
            file.write("\n\n")
