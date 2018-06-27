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

def symbol_all( nbit ):
    ret_list = []
    if nbit == 1:
        for a0 in range(0,2):
            ret_list.append(GFn.GFn([a0],1))
        return ret_list
    if nbit == 2:
        for a0 in range(0,2):
            for a1 in range(0,2):
                ret_list.append(GFn.GFn([a0,a1],2))
        return ret_list
    if nbit == 4:
        for a0 in range(0,2):
            for a1 in range(0,2):
                for a2 in range(0,2):
                    for a3 in range(0,2):
                        ret_list.append(GFn.GFn([a0,a1,a2,a3],4))
        return ret_list
    err_msg = "No symbol for nbit = ", str(nbit)
    raise ValueError(err_msg)

def find_characteristic( a ):
	i = 1
	product = a
	while 1:
		if product == a and i>1: return i
		product = product*a
		i = i+1

def assert_divisible( g_gfn, n ):
	g_int = [ int(g) for g in g_gfn ]
	r = gf2_remainder(np.array([1]+[0]*(n-1)+[1]), np.array(g_int))
	if np.count_nonzero(r):
		print("Remainder = ", r)
		raise ValueError(g_int, "is not divisible to x^n+1")

def alphas( m ):
	x_list = []
	i = 0
	while 1:
	# for i in range(0,2**m):
		content = np.zeros(2**m)
		content[i] = 1
		x = GFn.GFn( fit_gfn(content, m),m )
		if i>0 and x == x_list[0]: return x_list
		x_list.append(x)
		i += 1
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

def min_weight( g, n ):
	g = np.poly1d([s.toGF2(1) for s in g])
	# print("reducted g = ", g)
	symbols = []
	num_of_info = 2**(n-len(g))
	print("num_of_info = ", num_of_info)
	for value in range(1,num_of_info):
		bin_list = [int(ii) for ii in bin(value)[2:]]
		bin_list_aligned = [0]*(n-len(bin_list)) + bin_list
		gf_list = []
		for s in bin_list_aligned:
			if s: gf_list.append(GFn.GFn(np.array([1]),1))
			else: gf_list.append(GFn.GFn(np.array([0]),1))
		symbols.append( gf_list )

	min_weight = n
	for i, s in enumerate(symbols):
		y = np.polymul(s,g)
		# print("Info #", i, ":", s)
		# print("y = symbol * g(x) = ")
		# print(y)
		weight = sum([int(y2)!=0 for y2 in y])
		min_weight = min( weight, min_weight )
		# print("weight = ", weight)
		# print("------------")
	print("min_weight = ", min_weight)
	return min_weight

n = 5
m = find_m(n,2)
o = GFn.GFn(np.array([1]+[0]*(m-1)),m)
z = GFn.GFn(np.array([0]+[0]*(m-1)),m)
g = np.poly1d([o,o,o,o,o])
a = GFn.GFn(np.array([0,1]+[0]*(m-2)),m)
if find_characteristic(a) is not 2**m:
	print("char = ", find_characteristic(a))
	raise ValueError("Characteristic of alpha is not 2^m")
assert_divisible( g, n )

print("(n,m,q) = (", n, m, "2 )")
x_list = alphas(m)
for i, x in enumerate(x_list):
	print("alphas #", i, ": alpha ^", i, "= ", x)

print("g(x) = ")
print(g)
powers, roots, eqs = find_roots( x_list, g )
for i, rp in enumerate(zip(roots,powers)):
	root, power = rp
	print("roots #", i, "roots =", root, "= alpha ^", power)
for i, xy in enumerate(eqs):
	x, y = xy
	print("eqs #", i, "g(", x, ") = ", y)

min_weight(g,n)

def find_BCH( powers, n ):
	max_d0 = 2
	for b0 in powers:
		d0 = 2
		while 1:
			d0 = d0+1
			if (b0+d0-2)%n in powers:
				d0
				# print( (b0+d0-2)%n, " is in ", powers, "(b0,d0) = ", b0, d0 )
			else:
				d0 = d0-1
				break
		# print("For b0 = ", b0, ", d0 = ", d0)
		max_d0 = max(d0, max_d0)
	print("max d0 = ", max_d0)

def find_extBCH( powers, n ):
	max_d0 = 2
	for b0 in powers:
		for b0s in powers:
			if   b0s == b0: continue
			elif b0s >  b0: s = b0s-b0
			elif b0s <  b0: s = b0s-b0+n

			if gcd(n,s) > 1: continue
			# print("Pick s = ", s)
			for d0 in range(2,2*n):
				# d0 = d0+1
				if (b0+s*(d0-2))%n in powers:
					d0
					# print( (b0+s*(d0-2))%n, " is in ", powers, "(b0,d0) = ", b0, d0 )
				else:
					d0 = d0-1
					break
			# print("---")
			max_d0 = max(d0, max_d0)
	# 	print("For b0 = ", b0, ", d0 = ", d0)
	print("max d0 = ", max_d0)

def find_tzeng( powers, n ):
	large = n*n
	max_d0 = 2
	max_d0k0 = 2
	max_d0k0_param = 0
	for b0 in powers:
		# print("Pick b0 = ", b0)
		for b0s in powers:
			if    b0s == b0: continue
			else: s = (b0s-b0) % n

			if gcd(n,s) > 1: continue
			# print("\tPick s = ", s)
			for d0 in range(2,large):
				# print("\t(s,d0) = ", s, d0)
				if (b0+s*(d0-2))%n in powers:
					max_k0 = n
					for b0s2 in powers:
						if    b0s2 == b0+s*(d0-2): continue
						else: s2 = (b0s2-b0+s*(d0-2)) % n

						if max_k0 == 0: break
						if gcd(n,s2) >= d0: continue
						# print("\t\tPick s2 = ", s2)
						for k0 in range(0,max_k0):
							if (b0+s*(d0-2)+s2*k0)%n in powers:
								d0
								# print( "\t\t\t(b0,d0,k0) = ", b0, d0, k0 )
							else:
								k0 = k0-1
								break
						max_k0 = min( k0, max_k0 )
					# print("\t\tmax (d0,k0) = ", d0, max_k0)
					if d0 + max_k0 > max_d0k0:
						max_d0k0 = d0 + max_k0
						max_d0k0_param = (b0,s,d0,s2,k0,max_k0)
				else:
					d0 = d0-1
					break
			max_d0 = max(d0, max_d0)
		# print("---")
	# 	print("For b0 = ", b0, ", d0 = ", d0)
	# print("max d0 = ", max_d0)
	print("max d0k0 = ", max_d0k0, ", (b0,s,d0,s2,k0,max_k0) = ", max_d0k0_param)

find_BCH( powers, 2**m-1 )
find_extBCH( powers, 2**m-1 )
find_tzeng( powers, 2**m-1 )
# def BCH_bound( n, b, gen ):

