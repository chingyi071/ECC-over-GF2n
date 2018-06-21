import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from util import *

zero = np.zeros(2, dtype=int)

class BCJR:
    def __init__( self, p_mat ):
        self.b = 1
        self.n = p_mat.shape[1]
        self.k = p_mat.shape[1]-p_mat.shape[0]
        vex = [[zero]]
        edg = []
        for layer in range(0, self.n):
            vex_new = []
            edg_new = []
            for state_last in vex[-1]:
                for symbol in range(0,2**self.b):
                    # State calculation: s = last + symbol * H[i,:]
                    product = self.symbol_mul( symbol, p_mat[:,layer] )
                    add_v   = self.symbol_add( state_last, product )

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
                position = np.array([num_layer,-self.arr2int(v)])
                pos_dict[ state_name ] = position
        
        # Draw each edge in layer i
        for i, e_layer in enumerate(self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((array2str(v0,i+start), array2str(v1,i+start+1), {'weight':a}))

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
        vex_nonzero = [x for x in self.vex[-1] if not np.array_equal(x,zero)]
        self.edg.append([(zero,0,zero)])
        for i in range( len(self.edg)-2, -1, -1 ):
            vex_next = vex_in_edg( self.vex[i+1], self.edg[i+1])
            edg_new_list = []
            for edge in self.edg[i]:
                v0, a, v1 = edge
                if in_list( vex_next, v1 ):
                    edg_new_list.append(edge)
            self.edg[i] = edg_new_list

    def symbol_mul( self, a, b ):
        return np.remainder( a*b, 2**self.b )

    def symbol_add( self, a, b ):
        return np.remainder( a+b, 2**self.b )

    def arr2int( self, arr ):
        total = 0
        for digit in arr:
            total = total*2 + digit
        return total
