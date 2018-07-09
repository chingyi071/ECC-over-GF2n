from BCJR import BCJR, symbol_all
from util import *
import argparse
import GFn

parser = argparse.ArgumentParser(description="Flip a switch by setting a flag")
parser.add_argument('--input', default='test_data/test1.csv')
parser.add_argument('--plot_sections', type=int, nargs=2, default=[0,5])
parser.add_argument('--b', type=int, default=1)
args = parser.parse_args()

p_mat = read_mat(args.input)
zero_state = np.array( [GFn.GFn(0,args.b)]*int(p_mat.shape[1]) )
vex = [[zero_state]]
edg = []

symbol_np_arr = np.empty( [int(p_mat.shape[0]), int(p_mat.shape[1])], dtype=GFn.GFn)
for x in range(0,p_mat.shape[0]):
    for y in range(0,p_mat.shape[1]):
        symbol_np_arr[x][y] = GFn.GFn( p_mat[x][y], args.b )


for layer in range(0, int(p_mat.shape[0])):
    vex_new = []
    edg_new = []
    symbol_layer = symbol_np_arr[layer]
    for v_last in vex[-1]:
        for symbol in symbol_all(args.b):
            # State calculation: s = last + symbol * H[i,:]
            add_v   = symbol * symbol_layer + v_last

            # Create new edge and append it
            edge = ( v_last, symbol, add_v )
            edg_new.append(edge)

            # Add this vertice into vertex list if not exist
            if not in_list( vex_new, add_v ):
                vex_new.append(add_v)
    vex.append(vex_new)
    edg.append(edg_new)

bcjr1 = BCJR(p_mat, b=args.b, vex=vex, edg=edg)
bcjr1.remove_nonzero()
bcjr1.plot_sections(args.plot_sections)
