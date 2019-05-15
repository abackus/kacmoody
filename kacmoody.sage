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
    R = self.root_system()
    L = R.root_lattice()
    S = L.simple_roots()
    self.multiplicities = {} # Could imagine checking if it is defined first, so as
                             # not to waste time calculating things that have
                             # already been calculated.
    D = self.symmetrized_matrix().dense_matrix()

    def c(beta):
        dense = beta.dense_coefficient_list()
        f = common_divisors(dense)
        to_return = 0
        for e in f:
            if beta/e in self.multiplicities:
                to_return = to_return + (float(1)/float(e))*self.multiplicities[beta/e]
        return to_return

    def c_to_root(coeff):
        #IT IS 0 INDEXED FOR INFINITE, 1 INDEXED FOR FINITE
        #THIS IS A BUG IN SAGE
        #return sum([S[x+1]*coeff[x] for x in range(len(coeff))])#finite
        return sum([S[int(x)]*int(coeff[x]) for x in range(len(coeff))])#finite

    #Actually, it is SO much easier to induct on height...
    all_combs = []
    for l in [all_lists(len(S),x) for x in range(1,depth+1)]:
        all_combs = all_combs + l
    roots = [c_to_root(comb) for comb in all_combs]


    def two_rho(beta):
        return -2*sum(beta.dense_coefficient_list())


    for r in roots:
        rd = r.dense_coefficient_list()
        if sum(rd) == 1:
            self.multiplicities[r] = 1
        else:
            a_o = 0
            for bp in [c_to_root(l) for l in sub_parts(rd) if sum(l) != 0]:
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

    while true:
        for i in range(dim):
            for j in range(dim):
                if (i != j): #Don't replace 2
                    matrix[i][j] = matrix[i][j] - 1
                    if (matrix[j][i] == 0):
                        matrix[j][i] = -1
                    yield CartanMatrix(matrix)

def numerology(dim, count, dep):
    # Computes root multiplicities for dimension dim with depth dep, until count Kac-Moody algebras have been studied.
    with open('data.txt', 'w') as file:
        gen = CartanMatrixGenerator(dim)
        for i in range(count):
            file.write(next(gen).rm(depth=dep)) #TODO debug
