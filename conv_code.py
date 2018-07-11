import GFn
import util
import numpy as np
from math import log2
from BCJR import BCJR, symbol_all
import argparse

def parse_golden( csv_name, logq, m, n ):
	output_csv = util.read_csv(csv_name)[0]
	golden_list = []
	for o in output_csv:
		o_symbol = [ GFn.GFn(int(o_char),logq) for o_char in util.sep( int(o), n, 2**logq) ]
		golden_list.append(o_symbol)
	golden_list.extend( [[GFn.GFn(0,logq)]*2]*m )
	golden_snap = np.array(golden_list)
	for i, o in enumerate(np.swapaxes(golden_snap,0,1)):
		print("output[", i, "] = ", o[::-1], "\n", GFn.GFn_poly( o[::-1], logq) )
	print("---")
	return golden_snap, golden_snap.shape[0]

def parse_gens( csv_name, logq, k, n, m ):
	g_csv = util.read_csv(csv_name)
	if not np.array_equal(g_csv.shape, (k,n)):
		err_msg = "Expected generators is (" + str(k) + ',' + str(n) + '), "' + csv_name + '" is ' + str(g_csv.shape)
		# print("g_csv.shape = ", g_csv.shape)
		# print("Ideal is ", k, n)
		raise ValueError(err_msg)
	gens = np.empty( shape=(k,n), dtype=object )
	for i in range( g_csv.shape[0] ):
		for j in range( g_csv.shape[1] ):
			# gens[i][j] = GFn.GFn_poly(bin(g_csv[i][j])[2:], logq)
			gens[i][j] = GFn.GFn_poly(util.sep( g_csv[i][j], m, 2**logq), logq)
			print("generators[", i, "][", j, "] = \n", gens[i][j])
	weights = np.empty( shape=(n,m,k), dtype=object )
	for nn in range(n):
		for kk in range(k):
			coeffs = util.zero_padding_front(gens[kk][nn].c, m, zero=GFn.GFn(0,logq))
			for mm in range(m):
				weights[nn][mm][kk] = coeffs[m-mm-1]
	print("---")
	return gens, weights

def bit_diff( a_arr, b_arr ):
	num_diff = 0
	for a, b in zip(a_arr,b_arr):
		if not (a+b).iszero():
			num_diff += 1
	return num_diff

def main():

	parser = argparse.ArgumentParser(description="Flip a switch by setting a flag")
	parser.add_argument('--verbose', action='store_true')
	parser.add_argument('--out_seq', default="conv_csv/output.csv", help="File path of output sequence")
	parser.add_argument('--gen', default="conv_csv/g.csv", help="File path of generators")
	parser.add_argument('--k', default=1, type=int, help='Number of input')
	parser.add_argument('--n', default=2, type=int, help='Number of output')
	parser.add_argument('--m', default=3, type=int, help='Number of memories')
	parser.add_argument('--q', default=2, type=int, help='GF(q)')
	args = parser.parse_args()

	k = args.k # num of input
	n = args.n # num of output
	m = args.m
	q = args.q
	N = 3
	logq = int(log2(q))

	# Define dimension of memories
	mems = np.reshape( [GFn.GFn(0,logq)]*(m*k), (m,k) )

	# Define all the possible input, which is a k-tuple over GF(q) => q^k kinds
	poss_input = []
	for int_value in range(q**k):
		poss_input.append([ GFn.GFn(x,logq) for x in util.sep( int_value, k, q) ])

	# Read given generators and output sequence from csv
	golden_snap, output_len = parse_golden( args.out_seq, logq=logq, m=m, n=n)
	gens, weights = parse_gens( args.gen, logq=logq, k=k, n=n, m=m)

	# Initialize vertexs and edges of the graph
	zero_state = np.array( [GFn.GFn(0,logq)]*m*k )
	vex = [[(zero_state,0)]]
	edg = []

	for i in range( output_len ):

		vex_new = []
		edg_new = []

		for last_mems_flat, last_diff in vex[-1]:
			# print("last_mems_flat = ", last_mems_flat)
			last_mems = np.reshape( last_mems_flat, (m,k))


			for input_bit in poss_input:
			# for input_bit in symbol_all(k):
				# print("input_bit = ", input_bit, type(input_bit))
				input_snap = np.reshape( input_bit, (1,k) )

				# Shift memory a column
				mems = np.concatenate( (input_snap, last_mems[:-1]), axis=0)
				# print("mem = ",  mems.shape, "\n", mems)

				# Calculate every output
				outputs = np.empty(n, dtype=object)
				for j in range(n):
					output_sum = GFn.GFn(0,logq)
					for kk in range(k):
						for mm in range(m):
							output_sum += mems[mm][kk] * weights[j][mm][kk]
							# print("j(",j,"),m(",m,"),k(",k,") output_sum += ",  mems[mm][kk] * weights[j][mm][kk])
							# print("mems[mm][kk] = ", mems[mm][kk], type(mems[mm][kk]))
							# print("weights[j][mm][kk] = ", weights[j][mm][kk], type(weights[j][mm][kk]))
					outputs[j] = output_sum
					# outputs[i][j] = output_sum
					# print("outputs[", j, "] = ", outputs[j])
				
				# output_total = 0
				# for i, o in enumerate(outputs):
				# 	output_total = output_total*2 + int(o)
				# print("output_total = ", output_total, outputs)

				# print("mem flat = ", mems.flatten())
				# print("output/golden snap", outputs, golden_snap)
				num_diff = bit_diff(outputs, golden_snap[i])
				# print("diff = ", num_diff)
				mems_flat = mems.flatten()
				vex_new_names = [ v[0] for v in vex_new ]
				if not util.in_list( vex_new_names, mems_flat ):
					vex_new.append((mems_flat, last_diff+num_diff))
					edg_new.append((last_mems_flat, last_diff+num_diff, mems_flat))
				else:
					# Find target vertex to be replaced
					index_target = None
					for j, v in enumerate(vex_new):
						v_name, v_diff = v
						if np.array_equal(v_name, mems_flat):
							if index_target is not None: raise Exception
							index_target = j

					# Replace that vertex and edge if the new dist is smaller
					if last_diff+num_diff <= vex_new[index_target][1]:
						vex_new[index_target] = (mems_flat, last_diff+num_diff)
						edg_new = [e for e in edg_new if not np.array_equal(e[2],mems_flat)]
						edg_new.append((last_mems_flat, last_diff+num_diff, mems_flat))


				# vex_new.append(mems.flatten())
		vex.append( vex_new )
		edg.append( edg_new )

	# Filter only the first element (name) of the vertex as vertex names list
	v_names = []
	for v_layer in vex:
		v_names.append([])
		for v_name, v_diff in v_layer:
			v_names[-1].append(v_name)

	# for i, v_layer in enumerate(v_names):
	# 	print("Layer #", i)
	# 	for v_name in v_layer:
	# 		print("\tv_names = ", v_name)


	bcjr1 = BCJR( n=output_len, k=1, b=logq, vex=v_names, edg=edg, state_num=k*m)
	bcjr1.remove_disconnected()
	bcjr1.remove_nonzero()

	# for i, e_layer in enumerate(bcjr1.edg):
	# 	print("Edge layer #", i)
	# 	for e in e_layer:
	# 		print("e = ", e)

	input_predicted = np.empty( shape=(k,len(bcjr1.edg)), dtype=object )
	for i, e_layer in enumerate(reversed(bcjr1.edg)):
		if len(e_layer) is not 1:
			raise Exception
		e_inputs = e_layer[0][2][:k]
		for j, e_input in enumerate(e_inputs):
			input_predicted[j][i] = e_input
	
	input_predicted_poly = [GFn.GFn_poly(seq) for seq in input_predicted]

	for o_index, o_gens in enumerate(zip( np.swapaxes(golden_snap,0,1), np.swapaxes(gens,0,1))):
		o, gens = o_gens
		o_predicted_poly = GFn.GFn_poly( 0, logq )
		print("Output #", o_index)
		
		for i, g in enumerate(gens):
			print("input_predicted_poly[", i, "] = ", input_predicted[i])
			print("g[", i, "] = ", g.c)
			o_predicted_poly += input_predicted_poly[i] * g

		o_golden_poly = GFn.GFn_poly( o[::-1] )
		print("o_recv  = \n", o_predicted_poly)
		print("o_ideal = \n", o_golden_poly)
		print("Number of error = ", GFn.weight(o_golden_poly + o_predicted_poly))
		print("---")
	# print("golden_snap.shape = ", golden_snap.shape, golden_snap)

	bcjr1.plot_sections([0,output_len])

main()