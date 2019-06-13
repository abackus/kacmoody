from pytrees import AVLTree

CartanMatrix.multiplicities = dict()
CartanMatrix.cs = dict()

"""
AUXILIARY METHODS FOR AVL TREES
move to another file maybe?
"""
def as_list(self):
	"""
	Turns an AVLTree into a list, in the order that the entries appear
	in the AVLTree.
	"""
	l = []
	def asListHelper(node):
		if node is None:
			return []
		return asListHelper(node.left) + [node.val] + asListHelper(node.right)

	return asListHelper(self.root)
AVLTree.as_list = as_list

def insert_list(self, l):
	"""
	Adds a list into an AVLTree.
	"""
	for x in l:
		self.insert(x)
AVLTree.insert_list = insert_list

"""
OTHER AUXILIARY METHODS
"""
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

class root(tuple):
    """
    The class for roots of a Kac-Moody algebra.
    """
    #def __init__(self,list_form,multiplicity,c,subroots=None):
    def __init__(self, list_form):
        self.list_form    = list(list_form)
        self.vector_form  = vector(list_form)
        self.tuple_form   = tuple(list_form)
        self.coding       = prime_code(list_form)
    def __str__(self):
        return str(self.list_form)
    def __repr__(self):
        return str(self.list_form)
    def __eq__(left, right):
    	return left.coding == right.coding
    def __le__(left,right):
        return all(l <= r for l, r in zip(left.list_form, right.list_form))
    def __lt__(left,right):
        return all(l <= r for l, r in zip(left.list_form, right.list_form)) and not left.coding == right.coding
    def height(self):
        return sum(self.list_form)
    def set_multiplicity(self, mat, m):
    	mat.multiplicities[self.coding] = m
    def multiplicity(self, mat):
    	return mat.multiplicities[self.coding]
    def set_c(self, mat, c):
    	mat.cs[self.coding] = c
    def c(self, mat):
    	return mat.cs[self.coding] # REplace

class RootedCartanMatrix(CartanMatrix):
	"""
	The class for CartanMatrices decorated with more information.
	"""
	def __init__(self, height=30):
		super().__init__(self)
		self.height = height
		self.dim = self.nrows()
		self.simple_roots = [root(tuple([int(a == b) for a in range(dim)])) for b in range(dim)]
		self.roots = AVLTree()
		self.zero = root(tuple([0] * self.dim))
		for s in self.simple_roots:
			s.set_multiplicity(self, 1)
			s.set_c(self, 1)
		self.pingpong(self.simple_roots) # Generate real roots

	def B(self, a, b):
		"""
		Returns the bilinear product induced by self - does this work for nonsymmetric?
		"""
		return a * self * b

	def weyl(self, r, s):
		"""
		Acts the Weyl group representation of the simple root s on the root r r
		"""
		return root(r.vector_form - self.B(r.vector_form, s.vector_form) * s.vector_form)

	def pingpong(self, generators):
		"""
		Adjoin all the roots, that can be obtained
		by acting the Weyl group on the generators.
		"""
		self.roots.insert_list(generators)
		for g in generators:
			to_pingpong = [g]
			mult = g.multiplicity(self)
			while len(to_pingpong) != 0:
				next_root    = to_pingpong.pop()
				ponged       = [weyl(self, next_root, s) for s in self.simple_roots]
				qonged       = [p for p in ponged if p.height() <= self.height and not self.roots.search(p) and self.zero < p]
				to_pingpong += qonged
				for p in qonged:
					p.set_multiplicity(self, mult)
					self.roots.insert(p)


def exceptional(n, height=30):
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
    return RootedCartanMatrix(M, height)
