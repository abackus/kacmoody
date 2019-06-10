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
    
    #Josh : Yeah I think the function as described runs pretty much as quick as it can
    #I tested a variant; which fits into one-line:
    #    return [list(y) for y in cartesian_product([range(x+1) for x in c])][:-1]  #< I kinda like this one-liner...
    #but it is slower by a factor of ~3. Which is a shame.
    
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

def rm(self,height=10):
    """
    Generates a root multiplicity table for the fundamental chamber of a kac-moody algebra up to
    a certain height.

    Assume symmetric for now.
    """
    self.multiplicities = {tuple([0,0,0]) : 0}
    self.c = {tuple([0,0,0]) : 0}
    self.fund_mult = {} #multiplicities in the fundamental chamber
    mat = self.symmetrized_matrix()
    dim = mat.nrows()
    simples = all_lists(dim,1)

    def two_rho(beta):
        return -2*sum(beta)
    def B(a,b):
        return a*mat*b
    def is_fund(r):
        #checks if r is in the fundamental chamber
        for s in simples:
            if B(vector(r),vector(s)) > 0:
                return 0
        return 1
    def weyl(r,s):
        #for a root r and a simple s; ping pong r by s
        rv = vector(r); sv = vector(s)
        return tuple(rv - B(rv,sv)*sv)


    for root in simples:
        self.multiplicities[tuple(root)] = 1
        self.c[tuple(root)] = 1
    #Let's ping pong right now to see whats up with the homies
    to_pingpong = [root for root in simples]
    seen = [root for root in simples]
    while len(to_pingpong) != 0:
        #print(root,to_pingpong)
        el = to_pingpong.pop()
        ponged = [p for p in [weyl(el,s) for s in simples] if sum(p) <= height and p not in seen]
        for p in ponged:
            for q in p:
                tb = 0
                if q < 0:
                    tb = 1; break
            if tb == 1:
                break
            seen.append(p)
            self.multiplicities[tuple(p)] = 1
            self.c[tuple(p)] = 1
            to_pingpong.append(p)

    for h in range(2,height+1):
        for root in all_lists(dim,h):

            #check if we already know the mult
            if tuple(root) in self.multiplicities:
                #we already found this; this is a weyl reflection of something else
                #print(self.multiplicities)
                self.c[tuple(root)] = sum([self.multiplicities[tuple(vector(root)/n)]/n for n in divisors(gcd(root))])
                #self.multiplicities[tuple(root)] = m
                if is_fund(root):
                    self.fund_mult[tuple(root)] = self.multiplicities[tuple(root)]
                continue

            if B(vector(root),vector(root)) > 0:
                self.multiplicities[tuple(root)] = 0
                self.c[tuple(root)] = sum([self.multiplicities[tuple(vector(root)/n)]/n for n in divisors(gcd(root))])
                continue

            #if we don't already know the c
            RHS = 0 # Right-hand side of Peterson's recurrence formula
            root = vector(root)
            for subroot in sub_parts(list(root))[1:]:
                subroot = vector(subroot)
                RHS = RHS + B(subroot, root - subroot)*self.c[tuple(subroot)]*self.c[tuple(root - subroot)]
            LHSFactor = B(root, root) + two_rho(root)
            cRoot = sum([self.multiplicities[tuple(root/n)]/n for n in divisors(gcd(root))[1:]])

            if RHS == 0 and LHSFactor == 0:
                # Not a root
                self.multiplicities[tuple(root)] = 0
                self.c[tuple(root)] = cRoot
                continue

            self.multiplicities[tuple(root)] = RHS/LHSFactor - cRoot
            self.c[tuple(root)] = RHS/LHSFactor
            #print(self.multiplicities[tuple(root)])
            if is_fund(root):
                #print('wow')
                self.fund_mult[tuple(root)] = self.multiplicities[tuple(root)]

            #Now, we need to ping pong the root outwards by the Weyl group up to the height we care about
            to_pingpong = [root]
            seen = [root]
            while len(to_pingpong) != 0:
                #print(root,to_pingpong)
                el = to_pingpong.pop()
                ponged = [p for p in [weyl(el,s) for s in simples] if sum(p) <= height and p not in seen]
                for p in ponged:
                    for q in p:
                        tb = 0
                        if q < 0:
                            tb = 1; break
                    if tb == 1:
                        break
                    seen.append(p)
                    self.multiplicities[tuple(p)] = self.multiplicities[tuple(root)]
                    self.c[tuple(p)] = self.c[tuple(root)]
                    to_pingpong.append(p)

    G = [(x,-B(vector(x),vector(x)),self.fund_mult[x]) for x in sorted(self.fund_mult,key=lambda x : -B(vector(x),vector(x)))]
    #print(G)
    return G


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




def print_two(H):
    """
    prints a bunch of data, H, into two tables just like Kac
    """
    buff = 5

    #First, takes H and makes it into three list
    h0 = [str(h[0]) for h in H]; m0 = max([len(x) for x in h0] + [4])
    h1 = [str(h[1]/2) for h in H]; m1 = max([len(x) for x in h1] + [4])
    h2 = [str(h[2]) for h in H]; m2 = max([len(x) for x in h2] + [4])

    #Then, we can probably print a header
    #|m0 buff|m1 buff|m2 buff|buff|m0 buff|m1 buff|m2 buff|
    print('_'*(8+m0*2+m1*2+m2*2+buff*7))
    print('|root'+' '*(m0-4+buff) + '|norm' + ' '*(m1-4+buff) + '|mult' + ' '*(m2-4+buff) +'|'+' '*buff + '|root'+' '*(m0-4+buff) + '|norm' + ' '*(m1-4+buff) + '|mult' + ' '*(m2-4+buff)+'|')
    print('-'*(8+m0*2+m1*2+m2*2+buff*7))

    m = ceil(len(H)/2)
    for a in range(m):
        a1='|'+h0[a]+' '*(m0-len(h0[a])+buff)+'|'+h1[a]+' '*(m1-len(h1[a])+buff)+'|'+h2[a]+' '*(m2-len(h2[a])+buff)+'|'
        a2=' '*buff
        a3='|'+h0[a+m]+' '*(m0-len(h0[a+m])+buff)+'|'+h1[a+m]+' '*(m1-len(h1[a+m])+buff)+'|'+h2[a+m]+' '*(m2-len(h2[a+m])+buff)+'|'
        print(a1+a2+a3)
    print('-'*(8+m0*2+m1*2+m2*2+buff*7))


def exceptional(n):
    # Returns the exceptional Cartan matrix E_n, for n >= 5
    if n < 5:
        raise ValueError("Matrix E_n only makes sense for n >= 5")
    M = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(2)
            elif (i == j + 1) or (j == i + 1):
                row.append(-1)
            elif ((i == n - 1) and (j == 2)) or ((j == n - 1) and (i == 2)):
                row.append(-1)
            else:
                row.append(0)
        M.append(row)
    M[n-1][n-2] = 0
    M[n-2][n-1] = 0
    return CartanMatrix(M)



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
