import numpy as np
import GFn
from util import gf2_remainder, fit_gfn, sep
from fractions import gcd
import math
import argparse

def find_m( n, q ):
	m = 1
	qm = q
	while 1:
		if (qm-1) % n == 0: return m
		m += 1
		qm *= q

def get_min_weight( g, n, verbose=0 ):

	k = n-(g.order)
	nbit = g[0].nbit
	num_of_info = (2**nbit)**k
	print("n,g,k,num_of_info = ", n,g.order,k,num_of_info)

	if verbose:
		print("reducted g = ", g)
		print("num_of_info = ", num_of_info)
	
	# Create symbols from 1 to 2^k
	symbols = []
	for value in range(1,num_of_info):
		bin_list_aligned = sep( value, n, 2**nbit)
		# print("(", bin_list_aligned, bin_list_aligned2, ")")
		gf_list = [ GFn.GFn(s,nbit) for s in bin_list_aligned ]
		symbols.append( GFn.GFn_poly(gf_list) )

	# Find min weight in the product of g(x) and all symbols
	min_weight = n+1
	if not len(symbols) == (num_of_info-1):
		print("symbols = ", len(symbols))
		print("num_info-1 = ", num_of_info-1)
		raise Exception
	for i, s in enumerate(symbols):
		y = s*g
		weight = y.weight()
		# y = np.polymul(s,g)
		# weight = GFn.weight(y)
		if weight < min_weight:
			min_weight = weight
			min_codeword = s
		if verbose:
			print("Info #", i, ":")
			print("y = symbol * g(x) = ")
			print(y)
			print("weight = ", weight)
			print("------------")

	if verbose: print("min_weight = ", min_weight)
	return min_weight, min_codeword

def find_BCH( n, ext, g, verbose=False ):

	zero_ext, one_ext, alpha_ext = GFn.gen_zero_one_alpha_overGFq(2**ext)
	alphas_power = [alpha_ext.power(i) for i in range(2**ext-1)]
	roots_power, roots = GFn.finding_roots( g, ext=log_ext, alpha=alpha_ext )

	max_d0 = 2
	large = n*n

	# Select a root as alpha^b0
	for b0 in roots_power:
		# Find max d0 in root lists
		d0 = 2
		while (b0+d0-2) % (2**ext-1) in roots_power:
			d0 = d0+1
		max_d0 = max(d0-1, max_d0)
	return max_d0

def find_extBCH( n, ext, g, verbose=False ):

	zero_ext, one_ext, alpha_ext = GFn.gen_zero_one_alpha_overGFq(2**ext)
	alphas_power = [alpha_ext.power(i) for i in range(2**ext-1)]
	roots_power, roots = GFn.finding_roots( g, ext=log_ext, alpha=alpha_ext )

	max_d0 = 2
	large = n*n
	# Select a root as alpha^b0
	for b0 in roots_power:
		# Select another root as alpha^(b0+s)
		for b0s in [_b0s for _b0s in roots_power if _b0s is not b0]:
			s = (b0s-b0) % (2**ext-1)
			if gcd(n,s) > 1: continue
			if verbose: print("\tPick s = ", s)

			# Find max d0 in root lists
			d0 = 2
			while (b0+s*(d0-2)) % (2**ext-1) in roots_power:
				d0 = d0+1
			max_d0 = max(d0-1, max_d0)
	return max_d0

def find_tzeng( n, ext, g, verbose=False, check=True ):

	zero_ext, one_ext, alpha_ext = GFn.gen_zero_one_alpha_overGFq(2**ext)
	alphas_power = [alpha_ext.power(i) for i in range(2**ext-1)]
	roots_power, roots = GFn.finding_roots( g, ext=log_ext, alpha=alpha_ext )

	large = n*n
	max_d0k0 = 2
	max_d0k0_param = (0,0,0,0,0)

	# Select a root as alpha^b0
	for b0 in roots_power:
		if verbose: print("Pick b0 = ", b0)

		# Select another root as alpha^(b0+s)
		for b0s in [_b0s for _b0s in roots_power if _b0s is not b0]:

			s = (b0s-b0) % (2**ext-1)
			if gcd(n,s) > 1: continue
			if verbose: print("\tPick s = ", s)

			# Increase d0 and try to find max (d0+k0)
			d0 = 2
			while (b0+s*(d0-2)) % (2**ext-1) in roots_power:

				# Select another root as alpha^(b0+s+s2)
				if verbose: print("\t(b0,s,d0) = ", b0, s, d0)
				for b0s2 in [_b0s2 for _b0s2 in roots_power if (_b0s2 is not b0 and _b0s2 is not b0+s*(d0-2)) ]:

					s2 = (b0s2-b0) % (2**ext-1)
					if gcd(n,s2) >= d0: continue
					if verbose: print("\t\tPick s2 = ", s2)

					# Find max k0 in root lists and valid k0
					i2 = 0
					all_in_power = True
					while all_in_power:
						for i1 in range(0,d0-1):
							if (b0+s*i1+s2*i2) % (2**ext-1) not in roots_power:
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
	if verbose:
		print("max d0k0 = ", max_d0k0, ", (b0,s,d0,s2,k0) = ", max_d0k0_param)

	if check:
		b0, s, d0, s2, k0 = max_d0k0_param
		for i1 in range(0, d0-1):
			for i2 in range(0, k0+1):
				if (b0+s*i1+s2*i2) % (2**ext-1) not in roots_power:
					print("(b0, s, i1, s2, i2) = ", b0, s, i1, s2, i2)
					raise Exception

	return max_d0k0

def find_conjugate( value, base, ext ):
	conjugates = []
	for index in range(0, 2**ext):
		e = (value * (2**base)**index) % (2**ext-1)
		if e not in conjugates:
			conjugates.append(e)
		else:
			return conjugates


def find_generators( n, m, q=2, verbose=0 ):
	log_q = int(math.log2(q))
	log_ext = m * log_q
	alpha = GFn.GFn(2,log_ext)

	if verbose:
		print("Step 1: Find all cyclotonics group and its generator")
	used = []
	cyclotonics = []
	generator_base = []
	for i in range(0,(q**m)-1):
		if i not in used and int(alpha.power(i).power(n))==1:
			conjugate_power = find_conjugate( i, base=log_q, ext=log_ext )
			conjugate = []
			poly = GFn.GFn_poly(1,log_ext)
			# print("poly = \n", poly)
			for i in conjugate_power:
				# poly *= term = alpha^i
				term = GFn.GFn_poly([GFn.GFn(1,log_ext),alpha.power(i)])
				# poly = np.polymul(poly,term)
				poly = poly*term
				# print("poly = \n", poly)
			generator_base.append(poly)
			cyclotonics.append(conjugate_power)
			used.extend(conjugate_power)

			if verbose:
				cyc, gen = conjugate_power, poly
				print("#", i, "cyclotonics = ", cyc)
				print("generator = \n", gen)
				print("----")

	if verbose:
		print("Step 3: Map generator polynomial from GF(2^",log_ext,") to GF(2^", int(math.log2(q)), ")")

	gens_gfm = []
	for i, base in enumerate(generator_base):
		print("Base = \n", base)
		gen_gfm = base.map_to(log_ext)
		gens_gfm.append( gen_gfm )
		if verbose:
			print("Irreducible base #"+str(i))
			print("Base coeff from ascending order = ", base.c)
			print("Final coef from ascending order = ", gen_gfm.c, "\n", gen_gfm)
			print("---")

	if verbose:
		print("Step 4: Create every generator polynomial from those irreducible term")
	generators = []
	for selection in range(1,2**len(gens_gfm)-1):
		bin_list = [int(ii) for ii in bin(selection)[2:]]
		bin_list_aligned = [0]*(len(gens_gfm)-len(bin_list)) + bin_list
		g = GFn.GFn_poly(1,log_ext)
		for sel, base in zip( bin_list_aligned, gens_gfm ):
			if sel: g *= base
		generators.append(g)
		if verbose:
			print("Generator #", str(selection), ": ")
			print("bin_list = ", bin_list_aligned)
			print("g = ", g.c)
			print(g)
			print("---")

	return generators

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Flip a switch by setting a flag")
	parser.add_argument('--verbose', action='store_true')
	parser.add_argument('--gen')
	parser.add_argument('--n', type=int, default=5)
	parser.add_argument('--q', type=int, default=4)
	args = parser.parse_args()

	n = args.n
	q = args.q
	m = find_m(n,q)
	log_q = int(math.log2(q))
	log_ext = m * log_q

	alpha_ext = GFn.GFn(2,log_ext)
	if GFn.find_characteristic(alpha_ext) is not q**m:
		print("char = ", GFn.find_characteristic(alpha_ext))
		raise ValueError("Characteristic of alpha is not 2^m")

	if args.gen is not None:
		gens = [GFn.GFn_poly( args.gen, log_q)]
	else:
		gens = find_generators(n,m,q,verbose=args.verbose)

	zero_q, one_q, alpha_q = GFn.gen_zero_one_alpha_overGFq(q)
	# xn_1 = np.array([one_q]+[zero_q]*(n-1)+[one_q])
	xn_1 = (GFn.GFn_poly(1,log_q) << n) + GFn.GFn_poly(1,log_q)
	for gg in gens:			
		# r = GFn.gfn_array_modulo( xn_1, gg )
		r = xn_1 % gg
		if r.weight() is not 0:
			raise Exception
	
	print("(n,m,q) = (", n, m, q, ")")
	x_list = [alpha_ext.power(i) for i in range(q**m-1)]
	if args.verbose:
		for i, x in enumerate(x_list):
			print("alphas #", i, ": alpha ^", i, "= ", x)

	for g in gens:
		roots_power, roots = GFn.finding_roots( g, ext=log_ext, alpha=alpha_ext )
		if len(roots) is not g.order:
			print("g(x) = ")
			print(g)
			err_msg = "g(x) is a ", str(g.order), "-degree polynomial, but it has", str(len(roots)), "roots"
			raise ValueError(err_msg)

		if args.verbose:
			for i, rp in enumerate(zip(roots,roots_power)):
				root, power = rp
				print("roots #", i, "roots =", root, "= alpha ^", power)

		print("g(x) = \n", g)
		print("Roots = ", roots_power)
		BCH_bound    = find_BCH(    g=g, n=n, ext=log_ext, verbose=args.verbose )
		extBCH_bound = find_extBCH( g=g, n=n, ext=log_ext, verbose=args.verbose )
		tzeng_bound  = find_tzeng(  g=g, n=n, ext=log_ext, verbose=args.verbose )
		print("BCH bound    = ", BCH_bound)
		print("extBCH bound = ", extBCH_bound)
		print("tzeng bound  = ", tzeng_bound)
		# def BCH_bound( n, b, gen ):

		w, uw = get_min_weight( g, n, verbose=args.verbose )
		print("Min weight   = ", w)
		print("Min codeword = ", uw)
		print()