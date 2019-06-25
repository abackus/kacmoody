def sub_parts(c):
	"""
	For a list c = [1,2,0], return
				  [[0,0,0],
				   [1,0,0],
				   [0,1,0],
				   [1,1,0],
				   [0,2,0],
				   [1,2,0]]
	"""
	if len(c) == 1:
		return [[x] for x in range(c[0])]

	if c[0] == 0:
		return [[0] + y for y in sub_parts(c[1:])]

	d = sub_parts(c[1:])
	g = sub_parts([c[0]-1]+c[1:])+[[c[0]]+y for y in sub_parts(c[1:])]+[[c[0]-1]+c[1:]]
	return g+[c]

class root(tuple):
	"""
	The class for roots of a Kac-Moody algebra.
	"""
	def __init__(self, list_form):
		self.list_form	= list(list_form)
		self.vector_form  = vector(list_form)
		self.tuple_form   = tuple(list_form)
		self.height	   = sum(self.list_form)
	def __str__(self):
		return 'root'+str(self.list_form)
	def __repr__(self):
		return 'root'+str(self.list_form)
	def __le__(left,right):
		return all(l <= r for l, r in zip(left.list_form, right.list_form))
	def __lt__(left,right):
		return all(l <= r for l, r in zip(left.list_form, right.list_form)) and not left == right
	def __mul__(self,other):
		return self.vector_form * other
	def __rmul__(self,other):
		return other*self.vector_form
	def __sub__(self,other):
		return self.vector_form - other
	def __rsub__(self,other):
		return other - self.vector_form
	def __div__(self,other):
		return self.vector_form/other
	def __truediv__(self,other):
		return self.vector_form/other
	def divisors(self, j):
		# Return all roots beta such that m * beta = self, for m >= j
		return [(root(self/m), m) for m in divisors(gcd(self.list_form)) if m >= j]
	def is_positive(self):
		return all([i >= 0 for i in self.list_form]) and any([i > 0 for i in self.list_form])
	def base_root(self):
		# Return the root of least height in the span of self
		return root(tuple(self/divisors(gcd(self.list_form))[-1]))

class KacMoodyAlgebra():
	"""
	The class for CartanMatrices decorated with more information.
	"""
	def __init__(self, matrix, height=30):
		#Setup
		self.matrix = CartanMatrix(matrix)
		self.height = height
		self.dim = self.matrix.nrows()
		self.simple_roots = [root(tuple([int(a == b) for a in range(self.dim)])) for b in range(self.dim)]
		self.roots = set(self.simple_roots)
		self.multiplicities = dict()
		self.cs = dict()
		self.zero_root = root(tuple([0] * self.dim))
		self.cs[self.zero_root] = 0
		
		for s in self.simple_roots:
			self.multiplicities[s] = 1
			self.cs[s] = 1
		print('Done generating simple roots...')
		self.pingpong(self.simple_roots)
		print('Done generating real roots...')
		print('Found ' + str(len(self.roots)) + ' real roots.')
		
		possible_imag_basis = [root(list(-x)) for x in Cone(matrix).dual().Hilbert_basis()]
		self.imag_basis = [r for r in possible_imag_basis if r > self.zero_root]
		print('Done generating Hilbert basis of imaginary roots...')
		self.graded_ascent()
		print('Done!')
		
		
	def B(self, a, b):
		"""
		Returns the bilinear product induced by self - does this work for nonsymmetric?
		"""
		return a * self.matrix * b

	def weyl(self, r, s):
		"""
		Acts the Weyl group representation of the simple root s on the root r r
		"""
		return root(r - self.B(r, s) * s)

	def pingpong(self, generators):
		"""
		Adjoin all the roots, that can be obtained
		by acting the Weyl group on the generators, up to height
		"""
		for g in generators:
			to_pingpong = [g]
			mult = self.multiplicities[g]
			c	= self.cs[g]
			while len(to_pingpong) != 0:
				next_root	= to_pingpong.pop()
				ponged	   = [self.weyl(next_root, s) for s in self.simple_roots]
				qonged	   = [p for p in ponged if p.height <= self.height and self.zero_root <= p]
				for p in qonged:
					if p not in self.roots:
						to_pingpong.append(p)
						self.multiplicities[p] = mult
						self.cs[p] = c
						self.roots.add(p)

	def peterson(self, r, subroot):
		# Adds to the RHS of the Peterson formula.
		# Iterate this method over all subroots of the Peterson formula.
		coroot = root(r - subroot)
		#print('hip')
		if coroot not in self.cs:
			self.cs[coroot] = sum([self.multiplicities[d[0]]/d[1] for d in coroot.divisors(1) if d[0] in self.roots])
		#print('hup')
		if self.cs[coroot] > 0:
			if subroot not in self.cs:
				self.cs[subroot] = sum([self.multiplicities[d[0]]/d[1] for d in subroot.divisors(1) if d[0] in self.roots])
			#print('hooray!')
			return self.cs[subroot] * self.cs[coroot] * self.B(subroot, coroot)
		return 0

	def graded_ascent(self):
		#First, generate all imaginary roots
		coeff_lists = sub_parts([floor(self.height/r.height) for r in self.imag_basis])[1:]
		#Only need coeffs up to a certain point to bound by height
		fun_chamber = []
		for x in coeff_lists:
			next_root = sum([self.imag_basis[i] * x[i] for i in range(len(x))])
			if sum(next_root) <= self.height:
				next_root = root(list(next_root))
				if next_root not in fun_chamber:
					fun_chamber.append(next_root)
		fun_chamber = sorted(fun_chamber, key=lambda r: r.height)

		for r in fun_chamber:
			# Now, compute the multiplicity and c-value of r
			RHS = 0 # right-hand side of the Peterson formula
			subroots = [s for s in self.roots if s < r] # could this be sped up?
			already_summed = set()
			for s in subroots:
				# Either a root is real (spacelike), or imaginary (lightlike or timelike).
				# The set of imaginary roots is very small, so we do not optimize for it.
				# If we know the behavior of a root, we know the behavior
				# of all roots in its span.
				# First, if s is an imaginary root, we do the dumbest thing possible.
				if self.B(s, s) <= 0:
					RHS = RHS + self.peterson(r, s)
				elif s not in already_summed:
					s = s.base_root()
					for n in range(1, floor(r.height/s.height) + 1):
						ns = root(n * s)
						if n > 1 and ns in subroots:
							already_summed.add(s)
						RHS = RHS + self.peterson(r, ns)
			self.cs[r] = RHS/(self.B(r, r) - 2 * r.height)
			print(str(self.B(r, r) - 2 * r.height))
			self.multiplicities[r] = self.cs[r] - sum([self.multiplicities[d[0]]/d[1] for d in r.divisors(2)])
			self.roots.add(r)
			self.pingpong([r])
			print(r, r.height, self.B(r, r), self.multiplicities[r], self.cs[r])
		print("Generated all imaginary roots...")


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
	return M