import numpy as np
import GFn
from util import gf2_remainder, fit_gfn
from fractions import gcd

def find_m( n, q ):
	m = 1
	qm = q
	while 1:
		if (qm-1) % n == 0: return m
		m += 1
		qm *= q

def assert_divisible( g_int, n ):
	r = gf2_remainder(np.array([1]+[0]*(n-1)+[1]), np.array(g_int))
	if np.count_nonzero(r):
		print("Remainder = ", r)
		raise ValueError(g_int, "is not divisible to x^n+1")

def alphas( m, rng ):
	x_list = []
	for i in rng:
		content = np.zeros(2**m)
		content[i] = 1
		x = GFn.GFn( fit_gfn(content, m),m )
		x_list.append(x)
	return x_list

def find_roots( x_list, g ):
	powers = []
	roots  = []
	eqs = []
	for e, x in enumerate(x_list):
		y = np.polyval(g,x)
		eqs.append((x,y))
		if y.iszero():
			powers.append(e)
			roots.append(x)
	return powers, roots, eqs

def get_min_weight( g, n, verbose=0 ):

	# Generate g(x) with coeff in GF(2)
	g = np.poly1d([s.toGF2(1) for s in g])
	k = n-len(g)
	num_of_info = 2**k
	if verbose:
		print("reducted g = ", g)
		print("num_of_info = ", num_of_info)
	
	# Create symbols from 1 to 2^k
	symbols = []
	for value in range(1,num_of_info):
		bin_list = [int(ii) for ii in bin(value)[2:]]
		bin_list_aligned = [0]*(n-len(bin_list)) + bin_list
		gf_list = [ GFn.GFn(s,1) for s in bin_list_aligned ]
		symbols.append( gf_list )

	# Find min weight in the product of g(x) and all symbols
	min_weight = n
	for i, s in enumerate(symbols):
		y = np.polymul(s,g)
		weight = sum([int(yy)!=0 for yy in y])
		min_weight = min( weight, min_weight )
		if verbose:
			print("Info #", i, ":", s)
			print("y = symbol * g(x) = ")
			print(y)
			print("weight = ", weight)
			print("------------")

	if verbose: print("min_weight = ", min_weight)
	return min_weight

def find_BCH( powers, n, verbose=False ):
	max_d0 = 2
	large = n*n

	# Select a root as alpha^b0
	for b0 in powers:
		# Find max d0 in root lists
		d0 = 2
		while (b0+d0-2)%n in powers:
			d0 = d0+1
		max_d0 = max(d0-1, max_d0)
	return max_d0

def find_extBCH( powers, n, verbose=False ):
	max_d0 = 2
	large = n*n
	# Select a root as alpha^b0
	for b0 in powers:
		# Select another root as alpha^(b0+s)
		for b0s in [_b0s for _b0s in powers if _b0s is not b0]:

			s = (b0s-b0) % n
			if gcd(n,s) > 1: continue
			if verbose: print("\tPick s = ", s)

			# Find max d0 in root lists
			d0 = 2
			while (b0+s*(d0-2))%n in powers:
				d0 = d0+1
			max_d0 = max(d0-1, max_d0)
	return max_d0

def find_tzeng( powers, n, verbose=False, check=True ):
	large = n*n
	max_d0k0 = 2
	max_d0k0_param = (0,0,0,0,0)

	# Select a root as alpha^b0
	for b0 in powers:
		if verbose: print("Pick b0 = ", b0)

		# Select another root as alpha^(b0+s)
		for b0s in [_b0s for _b0s in powers if _b0s is not b0]:

			s = (b0s-b0) % n
			if gcd(n,s) > 1: continue
			if verbose: print("\tPick s = ", s)

			# Increase d0 and try to find max (d0+k0)
			d0 = 2
			while (b0+s*(d0-2))%n in powers:

				# Select another root as alpha^(b0+s+s2)
				if verbose: print("\t(b0,s,d0) = ", b0, s, d0)
				for b0s2 in [_b0s2 for _b0s2 in powers if (_b0s2 is not b0 and _b0s2 is not b0+s*(d0-2)) ]:

					s2 = (b0s2-b0) % n
					if gcd(n,s2) >= d0: continue
					if verbose: print("\t\tPick s2 = ", s2)

					# Find max k0 in root lists and valid k0
					i2 = 0
					all_in_power = True
					while all_in_power:
						for i1 in range(0,d0-1):
							if (b0+s*i1+s2*i2)%n not in powers:
								all_in_power = False
								break
							elif verbose:
								print("\t\t\t(b0,s,i1,s2,i2) = ", b0, s, str(i1)+'/'+str(d0-2), s2, i2)
						if all_in_power: i2 = i2+1

					# Update k0 with i2-1, because i2 < k0
					k0 = i2-1
					if verbose: print("\t\t\t(d0,k0) = ", d0, k0)

					# Update max d0+k0 and store its config
					if d0 + k0 >= max_d0k0:
						max_d0k0 = d0 + k0
						max_d0k0_param = (b0,s,d0,s2,k0)

				if verbose: print("\t\tmax (d0,k0) = ", d0, k0)

				d0 = d0+1

		if verbose: print("---")
	print("max d0k0 = ", max_d0k0, ", (b0,s,d0,s2,k0) = ", max_d0k0_param)

	if check:
		b0, s, d0, s2, k0 = max_d0k0_param
		for i1 in range(0, d0-1):
			for i2 in range(0, k0+1):
				if (b0+s*i1+s2*i2)%n not in powers:
					print("(b0, s, i1, s2, i2) = ", b0, s, i1, s2, i2)
					raise Exception

	return max_d0k0

def find_conjugate( base, n, q=2 ):
	index = 1
	conjugates = []
	while 1:
		e = (base * 2**index) % n
		if e not in conjugates:
			conjugates.append(e)
		else:
			return conjugates
		index += 1

def find_conjugates( n, m, q=2, verbose=0 ):
	conjugates = []
	generator_base = []
	used = []
	for i in range(0,n):
		if i not in used:
			conjugate_power = find_conjugate( i, n, q )
			conjugate = []
			poly = np.poly1d([GFn.GFn(1,m)])
			for i in conjugate_power:
				content = np.zeros(2**m)
				content[i] = 1
				x = GFn.GFn( fit_gfn(content, m),m )
				conjugate.append(x)
				term = np.poly1d([GFn.GFn(1,m), GFn.GFn(fit_gfn(np.array([0]*i+[1]),m),m)])
				poly = np.polymul(poly,term)
				if verbose:
					print("x = ", x)
					print("term = ")
					print(term)
					print("poly_r = ")
					print(poly)
					print("poly_l = ")
					print(poly)
					for p in poly:
						print(p)
			generator_base.append(poly)
			conjugates.append(conjugate)
			used.extend(conjugate_power)
			if verbose: print("----------------")

	generators = []
	for selection in range(1,2**len(generator_base)-1):
		bin_list = [int(ii) for ii in bin(selection)[2:]]
		bin_list_aligned = [0]*(len(generator_base)-len(bin_list)) + bin_list
		g = np.poly1d([GFn.GFn(1,m)])
		for sel, base in zip( bin_list_aligned, generator_base ):
			if sel: g = np.polymul( g, base )
		if verbose:
			print("bin_list = ", bin_list_aligned)
			print("g = ")
			print(g)
		generators.append(g)

	return generators

if __name__ == "__main__":
	n = 7
	m = find_m(n,2)
	gen = None
	z = GFn.GFn(0,m)
	o = GFn.GFn(1,m)

	a = GFn.GFn(2,m)
	if GFn.find_characteristic(a) is not 2**m:
		print("char = ", find_characteristic(a))
		raise ValueError("Characteristic of alpha is not 2^m")

	import argparse
	parser = argparse.ArgumentParser(description="Flip a switch by setting a flag")
	parser.add_argument('--verbose', action='store_true')

	args = parser.parse_args()

	if gen:
		g_int = [1,1,1,1,1,1,1]
		assert_divisible( g_int, n )
		g = GFn.intlist_to_gfpoly( g_int, m )
		gens = [g]
	else:
		gens = find_conjugates(n,m)
	
	print("(n,m,q) = (", n, m, "2 )")
	for g in gens:

		print("generator = ")
		print(g)
		x_list = alphas(m, range(2**m-1))
		if args.verbose:
			for i, x in enumerate(x_list):
				print("alphas #", i, ": alpha ^", i, "= ", x)

		if args.verbose:
			print("g(x) = ")
			print(g)
		powers, roots, eqs = find_roots( x_list, g )
		if len(roots) is not g.order:
			print("g(x) = ")
			print(g)
			print("len(g) = ", len(g))
			err_msg = "g(x) is a ", str(len(g)-1), "-degree polynomial, but it has", str(len(roots)), "roots"
			raise ValueError(err_msg)

		if args.verbose:
			for i, rp in enumerate(zip(roots,powers)):
				root, power = rp
				print("roots #", i, "roots =", root, "= alpha ^", power)
			for i, xy in enumerate(eqs):
				x, y = xy
				print("eqs #", i, "g(", x, ") = ", y)


		BCH_bound    = find_BCH( powers, 2**m-1, verbose=args.verbose )
		extBCH_bound = find_extBCH( powers, 2**m-1, verbose=args.verbose )
		tzeng_bound  = find_tzeng( powers, 2**m-1, verbose=args.verbose )
		print("BCH bound    = ", BCH_bound)
		print("extBCH bound = ", extBCH_bound)
		print("tzeng bound  = ", tzeng_bound)
		# def BCH_bound( n, b, gen ):

		w = get_min_weight( g, n, verbose=args.verbose )
		print("Min weight   = ", w)
		print()