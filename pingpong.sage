from pytrees import AVLTree

def prime_code(list_form):
	"""
	Encode a list of integers as an integer using a "Godel code"
	by prime numbers.
	"""
	P = Primes()
	c = 1
	for i in range(len(list_form)):
		c *= P.unrank(i)^list_form[i]
	return c

class root():
    """
    The class for roots of a Kac-Moody algebra. 
    """
    def __init__(self,list_form,multiplicity,c,subroots=None):
        self.multiplicity = multiplicity
        self.c            = c
        self.list_form    = list_form
        self.vector_form  = vector(list_form)
        self.tuple_form   = tuple(list_form)
        self.subroots     = subroots
        self.coding       = prime_code(list_form)
    def __str__(self):
        return "this root : " + str(self.list_form)
    def __repr__(self):
        return "this root : " + str(self.list_form)
    def __le__(left,right):
        return all(l <= r for l, r in zip(left.list_form, right.list_form))
    def __lt__(left,right):
        return all(l <= r for l, r in zip(left.list_form, right.list_form)) and not left.coding == right.coding
    def height(self):
        return sum(self.list_form)
    def subroots(self):
    	# Returns the set of subroots, if it is known.
    	if self.subroots == None:
    		raise UserWarning("Tried to find unknown subroots!")
    	else:
    		return self.subroots

class T(tuple):
    """
    tuples; but enriched with some <= structure
    """
    def __le__(left, right):
        return all(l <= r for l, r in zip(left, right))
    def height(self):
        return sum(list(self))

"""
Auxiliary functions for the Cartan Matrix
"""

def B(mat,a,b):
    """
    Returning the bilinear product induced by the Cartan Matrix
    """
    return a*mat*b

def weyl(mat,r,s):
    """
    Weyl action on the root r by the simple root s
    """
    rv = r.vector_form; sv = s.vector_form
    return root(rv - B(mat,rv, sv)*sv, r.multiplicity, r.c) #is it c?


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

def real_roots(cartan,height):
    """
    Returns the poset of real roots for a given cartan matrix, up to the height specified.
    """
    rroots  = AVLTree()
    dim     = cartan.nrows()
    simples = [root(tuple([int(a == b) for a in range(dim)]), 1, 1) for b in range(dim)]
    zero    = root(tuple([0]*dim), 0, 0)
    rroots = AVLTree.buildFromList(simples)
    
    to_pingpong = [x for x in simples]
    while len(to_pingpong) != 0:
        el = to_pingpong.pop()
        pngd = [weyl(cartan,el,s) for s in simples] 
        qngd = [p for p in pngd if p.height() <= height and not rroots.search(p) and zero < p]
        to_pingpong += qngd
        #[rroots.add(p) for p in qngd]
        for p in qngd:
            rroots.insert(p)
        print(rroots.countNodes(), ";", len(to_pingpong))
        if rroots.countNodes() > 100:
        	break

    return rroots
