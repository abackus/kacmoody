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
			self.pingpong(s)
		print('Done generating real roots...')
		
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
		return root(r - self.B(r, s) * s) # Does this still work? TODO

	def pingpong(self, g):
		"""
		Adjoin all the roots, that can be obtained
		by acting the Weyl group on g, up to height
		"""
		a = gcd(g.list_form)
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

	def lookup_coroot_c(self, r):
		# Look up c of a coroot, where we might not know the answer
		if r in self.cs:
			return self.cs[r]
		else:
			d = gcd(r.list_form)
			if root(r/d) in self.roots:
				self.cs[r] = 1/d
				return 1/d
			return 0


	def graded_ascent(self):
		#First, generate all imaginary roots
		coeff_lists = sub_parts([floor(self.height/r.height) for r in self.imag_basis])[1:]
		#Only need coeffs up to a certain point to bound by height
		self.fun_chamber = []
		for x in coeff_lists:
			next_root = sum([self.imag_basis[i] * x[i] for i in range(len(x))])
			if sum(next_root) <= self.height:
				next_root = root(list(next_root))
				if next_root not in self.fun_chamber:
					self.fun_chamber.append(next_root)
		self.fun_chamber = sorted(self.fun_chamber, key=lambda r: r.height)
		print("Fundamental chamber generated...")

		# FOR TIME TESTING
		real_root_counter = 0
		imag_root_counter = 0

		for r in self.fun_chamber:
			# Now, compute the multiplicity and c-value of r
			#if self.dim == 2 and self.is_symmetric: # TODO: Doesn't work if nonsymmetric!
				# This is a nice optimization for the 2D symmetric case, wherein we can assume that r[1] < r[2].
			#	opposite = root([r.list_form[1], r.list_form[0]])
			#	if opposite in self.roots:
			#		self.roots.add(r)
			#		self.cs[r] = self.cs[opposite]
			#		self.multiplicities[r] = self.multiplicities[opposite]
			#		self.pingpong(r)
			#		continue
			# TODO: Possible optimization -- what if B(r, r) = 0? Does this imply that m(r) = rank(matrix)?
			# This is NOT the case in general, so when is it?
			RHS = 0 # right-hand side of the Peterson formula
			subroots = [s for s in self.roots if s < r] # could this be sped up?
			for s in subroots:
				coroot = root(r - s)
				norm = self.B(s, s)

				if norm > 0:
					# Real subroots
					pairing = self.B(s, r)
					current_sum = 0
					for n in range(1, floor(r.height/s.height) + 1):
						n_coroot = root(r - n * s)
						if not self.zero_root < n_coroot:
							break
						two_factor = n_coroot.height - (n * s.height)
						if two_factor < 0:
							continue
						current_sum = current_sum + (pairing - n * norm) * self.lookup_coroot_c(n_coroot) * (2 - (two_factor == 0))
						real_root_counter += 1
					RHS = RHS + self.multiplicities[s] * current_sum
				elif self.lookup_coroot_c(coroot) > 0:
					# Imaginary subroots
					two_factor = coroot.height - s.height
					if two_factor < 0:
						continue
					RHS = RHS +	self.cs[s] * self.cs[coroot] * self.B(s, coroot) * (2 - (two_factor == 0))
					imag_root_counter += 1
			
			self.cs[r] = RHS/(self.B(r, r) - 2 * r.height) # Does this still work??
			self.multiplicities[r] = self.cs[r] - sum([self.multiplicities[d[0]]/d[1] for d in r.divisors(2)])
			self.roots.add(r)
			self.pingpong(r)
		print(real_root_counter, imag_root_counter)
		print("Generated all imaginary roots...")

	def print_to_file(self, filename):
		# Prints fundamental chamber to filename.tsv, and
		# all roots to filename_raw.tsv
		with open(filename + '.tsv', 'w') as f:
			for r in self.fun_chamber:
				f.write(str(r) + '	' + str(r.height) + '	' + str(self.B(r, r)) + '	' + str(self.multiplicities[r]) + '	' + str(self.cs[r]) + '\n')
		#with open(filename + '_raw.tsv', 'w') as f:
		#	for r in self.roots:
		#		f.write(str(r) + '	' + str(r.height) + '	' + str(self.B(r, r)) + '	' + str(self.multiplicities[r]) + '	' + str(self.cs[r]) + '\n')


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