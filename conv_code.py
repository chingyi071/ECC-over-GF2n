import GFn
import util
import numpy as np
from math import log2
from BCJR import BCJR, symbol_all

def parse_golden( csv_name, logq, m ):
	output_csv = util.read_csv("output.csv")[0]
	golden_list = []
	for o in output_csv:
		o_symbol = [ GFn.GFn(int(o_char),logq) for o_char in util.sep( int(o), 2, 2) ]
		golden_list.append(o_symbol)
	golden_list.extend( [[GFn.GFn(0,logq)]*2]*m )
	golden_snap = np.array(golden_list)
	print("golden_snap.shape = ", golden_snap.shape, golden_snap)
	return golden_snap, golden_snap.shape[0]

def parse_gens( csv_name, logq, k, n, m ):
	g_csv = util.read_csv("g.csv")
	gs = np.empty( shape=(k,n), dtype=object )
	for i in range( g_csv.shape[0] ):
		for j in range( g_csv.shape[1] ):
			gs[i][j] = GFn.GFn_poly(bin(g_csv[i][j])[2:], logq)
	print("gs = ", gs)

	w = np.empty( shape=(n,m,k), dtype=object )
	for nn in range(n):
		for kk in range(k):
			for mm in range(m):
				w[nn][mm][kk] = gs[kk][nn].c[m-mm-1]
	return gs, w

def main():
	k = 1 # num of input
	n = 2 # num of output
	m = 3
	q = 2
	N = 3
	logq = int(log2(q))

	golden_snap, input_len = parse_golden("output.csv", logq=logq, m=m)
	gs, w = parse_gens("g.csv", logq=logq, k=k, n=n, m=m)
	mems = np.reshape( [GFn.GFn(0,logq)]*(m*k), (m,k) )

	zero_state = np.array( [GFn.GFn(0,logq)]*m )
	vex = [[(zero_state,0)]]
	edg = []
	for i in range( input_len ):

		vex_new = []
		edg_new = []
		for last_mems_flat, last_diff in vex[-1]:
			# print("last_mems_flat = ", last_mems_flat)
			last_mems = np.reshape( last_mems_flat, (m,k))

			# for input_bit in [GFn.GFn(0,logq), GFn.GFn(1,logq)]:
			for input_bit in symbol_all(1):
				input_snap = np.reshape( input_bit, (1,k) )
				mems = np.concatenate( (input_snap, last_mems[:-1]), axis=0)
				# print("mem = ",  mems.shape, "\n", mems)

				outputs = np.empty(n, dtype=object)
				for j in range(n):
					output_sum = GFn.GFn(0,logq)
					for kk in range(k):
						for mm in range(m):
							output_sum += mems[mm][kk] * w[j][mm][kk]
							# print("j(",j,"),m(",m,"),k(",k,") output_sum += ",  mems[mm][kk] * w[j][mm][kk])
							# print("mems[mm][kk] = ", mems[mm][kk], type(mems[mm][kk]))
							# print("w[j][mm][kk] = ", w[j][mm][kk], type(w[j][mm][kk]))
					outputs[j] = output_sum
					# outputs[i][j] = output_sum
					# print("outputs[", j, "] = ", outputs[j])
				
				# output_total = 0
				# for i, o in enumerate(outputs):
				# 	output_total = output_total*2 + int(o)
				# print("output_total = ", output_total, outputs)

				# print("mem flat = ", mems.flatten())
				# print("output/golden snap", outputs, golden_snap)
				num_diff = 0
				# print("outputs / golden_snap", i, " = ", outputs.shape, golden_snap[i].shape)
				for a, b in zip(outputs, golden_snap[i]):
					if not a == b: num_diff += 1
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

	v_names = []
	for v_layer in vex:
		v_names.append([])
		for v_name, v_diff in v_layer:
			v_names[-1].append(v_name)
	# print("v_names = ", v_names)


	bcjr1 = BCJR( n=input_len, k=1, b=1, vex=v_names, edg=edg, state_num=3)
	bcjr1.remove_disconnected()
	bcjr1.remove_nonzero()
	bcjr1.plot_sections([0,input_len])

	input_predicted = []
	for i, e_layer in enumerate(bcjr1.edg):
		# print("Edge layer #", i)
		for e in e_layer:
			# print("e = ", e)
			# print("e = ", e, ", input = ", e[2][0])
			input_predicted = input_predicted + [e[2][0]]
	input_predicted_poly = GFn.GFn_poly( input_predicted[::-1] )
	for o_index, o_gens in enumerate(zip( np.swapaxes( golden_snap, 0, 1 ), np.swapaxes(gs,0,1))):
		o, gens = o_gens
		o_predicted_poly = GFn.GFn_poly( 0, logq )
		print("Output #", o_index)
		print("input_predicted_poly = \n", input_predicted_poly)
		
		for g in gens:
			print("g = \n", g)
			o_predicted_poly += input_predicted_poly * g

		o_golden_poly = GFn.GFn_poly( o[::-1] )
		print("o_recv  = \n", o_predicted_poly)
		print("o_ideal = \n", o_golden_poly)
		print("Number of error = ", (o_golden_poly + o_predicted_poly).weight())
		print("---")
	# print("golden_snap.shape = ", golden_snap.shape, golden_snap)

main()