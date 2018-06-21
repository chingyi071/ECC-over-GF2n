import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from util import *

def symbol_all( nbit ):
    return [sbl(i,nbit) for i in range(0, 2**nbit)]

class sbl:
    def __init__( self, value, nbit ):
        self.nbit = nbit
        self.value = value

    def __add__( self, a ):
        return sbl( np.remainder( a.value + self.value, 2 ), self.nbit)

    def __mul__( self, a ):

        if type(a) is type(np.array([])):
            product = np.empty_like(a)
            for x in range(0,product.shape[0]):
                product[x] = self*a[x]
            return product

        return sbl( a.value * self.value, self.nbit)

    def __eq__( self, a ):
        if not self.nbit == a.nbit:
            raise Exception()
            return False
        if not self.value == a.value:
            return False
        return True

    def __repr__( self ):
        return str(self.value)

    def asint(self):
        total = 0
        for digit in self.value:
            total = total*2 + digit
        return total

    def __int__(self):
        return int(self.value)

    def __str__(self):
        return int2bistr(int(self.value), self.nbit)

zero = np.array([sbl(0,2), sbl(0,2)])

def sbl_arr2int( v ):
    total = 0
    for digit in v:
        # print("digit = ", digit)
        total = total*2 + int(digit)
    return total

class BCJR:
    def __init__( self, p_mat, b=1 ):
        self.b = b
        self.n = p_mat.shape[1]
        self.k = p_mat.shape[1]-p_mat.shape[0]
        self.h = p_mat.shape[0]

        symbol_np_arr = np.empty_like( p_mat, dtype=sbl)
        for x in range(0,symbol_np_arr.shape[0]):
            for y in range(0,symbol_np_arr.shape[1]):
                value = p_mat[x][y]
                symbol_np_arr[x][y] = sbl( value, self.b )

        vex = [[sbl(0,self.k) * zero]]
        edg = []
        for layer in range(0, self.n):
            vex_new = []
            edg_new = []
            for state_last in vex[-1]:
                for symbol in symbol_all(1):
                    # State calculation: s = last + symbol * H[i,:]
                    product = symbol * symbol_np_arr[:,layer]
                    add_v   = product + state_last

                    # Create new edge and append it
                    edge = (state_last,symbol,add_v)
                    edg_new.append(edge)

                    # Add this vertice into vertex list if not exist
                    if not in_list( vex_new, add_v ):
                        vex_new.append(add_v)
            vex.append(vex_new)
            edg.append(edg_new)
        self.vex = vex
        self.edg = edg
        self._remove_nonzero()

    def plot_sections( self, start, end ):
        bit_per_state = self.n - self.k
        edges = []
        pos_dict = {}

        # Draw each vertex in layer i
        for num_layer, v_layer in zip(range(start,end+1), self.vex[start:end+1]):
            for v in v_layer:
                state_name = str(num_layer) + '_' + arr2bistr(v,bit_per_state)
                position = np.array([num_layer,-sbl_arr2int(v)])
                pos_dict[ state_name ] = position
        
        # Draw each edge in layer i
        for i, e_layer in enumerate(self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((array2str(v0,i+start), array2str(v1,i+start+1), {'weight':str(a)}))

        G = nx.Graph()
        G.add_edges_from(edges)

        pos = nx.spring_layout(G) # positions for all nodes
        nx.draw_networkx_nodes(G,pos_dict,node_size=700,)
        nx.draw_networkx_labels(G,pos_dict,font_size=12,font_family='sans-serif')
        nx.draw_networkx_edges(G,pos_dict,edgelist=edges, width=6)
        labels = nx.get_edge_attributes(G,'weight')
        nx.draw_networkx_edge_labels(G,pos_dict,edge_labels=labels)
        plt.show()

    def _remove_nonzero( self ):
        zero_state = np.array( [sbl(0,self.b)] * self.h )
        self.edg.append([(zero_state,0,zero_state)])
        for i in range( len(self.edg)-2, -1, -1 ):
            vex_next = vex_in_edg( self.vex[i+1], self.edg[i+1])
            edg_new_list = []
            for edge in self.edg[i]:
                v0, a, v1 = edge
                if in_list( vex_next, v1 ):
                    edg_new_list.append(edge)
            self.edg[i] = edg_new_list
