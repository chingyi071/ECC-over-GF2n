import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from util import *
import GFn

def symbol_all( nbit, b=1 ):
    ret_list = []
    if nbit == 1:
        for a0 in range(0,2**b):
            ret_list.append(GFn.GFn([a0],1))
        return ret_list
    if nbit == 2:
        for a0 in range(0,2):
            for a1 in range(0,2**b):
                ret_list.append(GFn.GFn([a0,a1],2))
        return ret_list
    err_msg = "No symbol for nbit = ", str(nbit)
    raise ValueError(err_msg)

class BCJR:
    def __init__( self, n, k, vex, edg, state_num, b=1 ):
        self.b = b
        self.n = n
        self.k = k
        self.h = self.n - self.k

        zero_state = np.array( [GFn.GFn(0,self.b)]*state_num )
        self.zero_state = zero_state
        self.vex = vex
        self.edg = edg

    def plot_sections( self, position ):
        start, end = position
        bit_per_state = self.n - self.k
        edges = []
        pos_dict = {}

        # Draw each vertex in layer i
        for num_layer, v_layer in zip(range(start,end+1), self.vex[start:end+1]):
            for v in v_layer:
                v_name = v2str(v,num_layer)
                position = np.array([ num_layer, -arr2int(v) ])
                pos_dict[ v_name ] = position
        
        # Draw each edge in layer i
        for num_layer, e_layer in zip(range(start,end), self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((v2str(v0,num_layer), v2str(v1,num_layer+1), {'weight':int(a)}))

        G = nx.Graph()
        G.add_edges_from(edges)

        pos = nx.spring_layout(G) # positions for all nodes
        nx.draw_networkx_nodes(G,pos_dict,node_size=800)
        nx.draw_networkx_labels(G,pos_dict,font_size=12,font_family='sans-serif')
        nx.draw_networkx_edges(G,pos_dict,edgelist=edges, width=6)
        labels = nx.get_edge_attributes(G,'weight')
        nx.draw_networkx_edge_labels(G,pos_dict,edge_labels=labels)
        plt.show()

    def remove_nonzero( self ):
        self.edg.append([(self.zero_state,0,self.zero_state)])

        for i, e_layer in reversed( list( enumerate(self.edg[:-1]))):
            vex_next = vex_connected( self.vex[i+1], self.edg[i+1])
            edg_new_list = []
            for edge in e_layer:
                _, _, v1 = edge
                if in_list( vex_next, v1 ):
                    edg_new_list.append(edge)
            self.edg[i] = edg_new_list

    def remove_disconnected( self ):
        for i, e_layer in reversed( list( enumerate(self.edg[:-1]))):
            vex_next = vex_connected( self.vex[i+1], self.edg[i+1])
            edg_new_list = []
            for edge in e_layer:
                _, _, v1 = edge
                if in_list( vex_next, v1 ):
                    edg_new_list.append(edge)
            self.edg[i] = edg_new_list