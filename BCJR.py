import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from util import *

def prmt_ply( nbit ):
    if nbit==1: return np.array([0,1])
    if nbit==2: return np.array([1,1,1])
    d
def symbol_all( nbit ):
    ret_list = []
    if nbit == 1:
        for a0 in range(0,2):
            ret_list.append(sbl([a0],1))
        return ret_list
    if nbit == 2:
        for a0 in range(0,2):
            for a1 in range(0,2):
                ret_list.append(sbl([a0,a1],2))
        return ret_list
    raise ValueError()
    return [sbl(i,nbit) for i in range(0, 2**nbit)]

def np_remainder( a, b ):
    # print('a = ', a)
    # print('b = ', b)
    if a.size < b.size:
        # print('r = ', a[:b.size-1])
        return a[:b.size-1]
    while 1:
        if np.argwhere(a==1).size == 0:
            # print('r = ', a[:b.size-1])
            return a[:b.size-1]
        msb = np.max(np.argwhere(a==1),axis=0)[0]
        if msb < len(b)-1:
            # print('r = ', a[:b.size-1])
            return a[:b.size-1]
        remainder = np.zeros( msb - b.size + 1 )
        remainder = np.append( remainder, b )
        remainder = np.append( remainder, np.zeros(a.size - msb - 1))
        result = np.empty_like(a)
        for i, x in enumerate(result):
            result[i] = np.remainder( a[i]+remainder[i], 2)
        a = result

class sbl:
    def __init__( self, value, nbit ):
        self.nbit = nbit
        # print('value = ', value)
        if len(np.shape(value)) is not 1:
            raise ValueError()
        self.value = np.array(value).flatten().astype(int)

    def __add__( self, a ):
        if self.value.shape[0] is not a.value.shape[0]:
            print('self = ', self)
            print('a = ', a)
            print('self.value.shape = ', self.value.shape)
            print('a.value.shape = ', a.value.shape)
            raise ValueError()
        result_value = np.remainder(a.value + self.value,2)
        result_value = np_remainder(result_value.astype(int),prmt_ply(self.nbit))
        result = sbl( result_value, self.nbit)
        print(a, '+', self, '=', result)
        print(a.value, '+', self.value, '=', result.value)
        return result

    def __mul__( self, a ):

        # if type(a) is type(np.array([])):
        #     product = np.empty_like(a)
        #     for x in range(0,product.shape[0]):
        #         product[x] = self*a[x]
        #     return product
        if len(self.value.shape) > 1:
            raise ValueError
        if type(a) is type(np.array([])):
            flat_input = np.reshape( a, a.size )
            flat_output = np.zeros( a.size, dtype=object )
            for i in range(0, len(flat_input)):
                flat_output[i] = self*flat_input[i]
            return np.reshape(flat_output,a.shape)

        if type(a) is not type(sbl([0],1)):
            raise Exception
        product = np.zeros_like(a.value)
        # print('a.value = ', a.value)
        for i in range( 0, self.nbit ):
            z = np.zeros(i)
            # print('z1 = ', z, z.shape)
            z = np.append( z, self.value[i]*a.value)
            z = np.append( z, np.zeros(self.nbit-i))
            # print('z2 = ', z, z.shape)
            z1 = z[:self.nbit]
            z2 = np_remainder(z.astype(int),prmt_ply(self.nbit))
            # print('z1,z2 = ', z1, z2)
            z = z2
            if z.shape[0] is not a.nbit:
                print('mul (', self.value, ', ', a.value, ')')
                print('a.type = ', type(a), np.shape(a))
                print('z = ', z, z.shape)
                print('expected nbit = ', a.nbit)
                print('self.value = ', self.value)
                print('self.nbit = ', self.nbit)
                raise Exception()
            z = z.flatten().astype(int)
            # print('z = ', z, z.shape)
            # print('product = ', product, product.shape)
            product += z
        # print('product = ', product, product.shape)
        # print('z = ', z, z.shape)
        if product.shape[0] is not a.nbit:
            raise Exception()
        result = sbl( np.remainder(product, np.ones_like(a.value)*2), a.nbit)
        # print(a, '(', str(a), ")*", str(self), "=", str(result))
        return result
        # return sbl( product, a.nbit)
        # return sbl( a.value * value, self.nbit)

    def __eq__( self, a ):
        if not self.nbit == a.nbit:
            raise Exception()
            return False
        for index in range(0, self.nbit):
            if not self.value[index] == a.value[index]:
                return False
        return True

    def __repr__( self ):
        return 'sbl(' + str(self) + ')'

    def asint(self):
        total = 0
        for i, digit in enumerate(reversed(self.value)):
            total *= 2
            total += int(digit)
        return total

    def __int__(self):
        total = 0
        for i, digit in enumerate(reversed(self.value)):
            total *= 2
            total += int(digit)
        return total

    def __str__(self):
        out_str = ""
        for i, x in enumerate(self.value):
            out_str = str(x) + out_str
        # print('str -------------- ', self.value, out_str)
        out_str = str(int(self))
        return out_str

zero = np.array([sbl([0],1), sbl([0],1)])

def arr2int( v ):
    total = 0
    # print('v = ', v, v.shape)
    for e, digit in enumerate(v):
        # print("digit = ", digit)
        total = total*(2**digit.nbit) + int(digit)
    return total

class BCJR:
    def __init__( self, p_mat, b=1 ):
        self.b = b
        self.n = int(p_mat.shape[0])
        self.k = int(p_mat.shape[0])-p_mat.shape[1]
        self.h = p_mat.shape[1]

        # print('p_mat.size = ', p_mat.shape)
        # p_mat = np.reshape( p_mat, (int(p_mat.shape[0]/self.b), int(p_mat.shape[1]), self.b))
        # symbol_np_arr = np.empty_like( p_mat, dtype=sbl)
        symbol_np_arr = np.empty( [int(p_mat.shape[0]), int(p_mat.shape[1])], dtype=sbl)
        # for x in range(0,symbol_np_arr.shape[0]):
        #     for y in range(0,symbol_np_arr.shape[1]):
        #         value = p_mat[x][y]
        #         print('value = ', value)
        #         symbol_np_arr[x][y] = sbl( value, self.b )
        print('symbol_np_arr.size = ', symbol_np_arr.shape)
        print('p_mat.size = ', p_mat.shape)
        for i,x in enumerate(p_mat):
            print('p_mat[', i, '] = ', x)
        for x in range(0,p_mat.shape[0]):
            for y in range(0,p_mat.shape[1]):
                value = []
                for b in range(0,p_mat.shape[2]):
                    value.append(p_mat[x][y][b])
                print('value = ', value)
                symbol_np_arr[x][y] = sbl( value, self.b )
        zero_state = np.reshape(np.array([sbl([0]*self.b,self.b)]*int(self.h/self.b)), (int(self.h/self.b),))
        vex = [[zero_state]]
        edg = []
        # print('vex[0][0].value.shape = ', vex[0][0].value.shape, type(vex[0]))
        # raise ValueError()
        # if vex[0][0].value.shape[0] is not 2 and vex[0][0].value.shape[1] is not 1:
        #     print('vex[0][0].value.shape = ', vex[0][0].value.shape)
        #     raise ValueError()
        for layer in range(0, self.n):
            vex_new = []
            edg_new = []
            for state_last in vex[-1]:
                for symbol in symbol_all(self.b):
                    # State calculation: s = last + symbol * H[i,:]
                    product = symbol * symbol_np_arr[layer]
                    add_v   = product + state_last

                    # Create new edge and append it
                    edge = (state_last,symbol,add_v)
                    print('symbol = ', symbol, type(symbol))
                    print('symbol_np_arr[layer] = ', symbol_np_arr[layer], type(symbol_np_arr[layer]))
                    print('state_last = ', state_last, type(state_last))
                    print('add_v = ', add_v, type(add_v))
                    # print('edge = ', edge, type(edge))
                    # print('edge[0] = ', edge[0], type(edge[0]))
                    edg_new.append(edge)

                    # Add this vertice into vertex list if not exist
                    if not in_list( vex_new, add_v ):
                        vex_new.append(add_v)
                    # print('edg_new = ', edg_new, type(edg_new))
                    # print('edg_new[0] = ', edg_new[0], type(edg_new[0]))
                # raise ValueError
            vex.append(vex_new)
            edg.append(edg_new)
        self.vex = vex
        for i, v_row in enumerate(vex):
            for v in v_row:
                print("#v", i, " v = ", v)
        for i, e_row in enumerate(edg):
            for e in e_row:
                print("#e", i, " e = ", e)
        self.edg = edg
        self._remove_nonzero()
        # raise ValueError

    def plot_sections( self, start, end ):
        bit_per_state = self.n - self.k
        edges = []
        pos_dict = {}

        # Draw each vertex in layer i
        for num_layer, v_layer in zip(range(start,end+1), self.vex[start:end+1]):
            # print('v_layer(', len(v_layer), type(v_layer[0]), ') = ', v_layer)
            for v in v_layer:
                # state_name = str(num_layer) + '_' + arr2bistr(v,bit_per_state)
                state_name = array2str(v,num_layer)
                position = np.array([num_layer,-arr2int(v)])
                pos_dict[ state_name ] = position
        print("pos_dict = ", pos_dict)
        
        # Draw each edge in layer i
        # for i, e_row in enumerate(self.edg):
        #     print("#", i, " e = ", e_row)
        for i, e_layer in enumerate(self.edg[start:end]):
            for e in e_layer:
                v0, a, v1 = e
                edges.append((array2str(v0,i+start), array2str(v1,i+start+1), {'weight':int(a)}))
        print('edges = ', edges)

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
        vex = np.array([sbl([0]*self.b,self.b)]*self.h)
        zero_state = (vex,sbl([0],1),vex)
        # for i, v_row in enumerate(self.vex):
        #     print("#v", i, " v = ", v_row)        
        # for i, e_row in enumerate(self.edg):
        #     print("#e", i, " e = ", e_row)        
        # print('zero_state = ', zero_state)
        self.edg.append([(vex,sbl([0],1),vex)])
        # print('edg = ', self.edg, type(self.edg))
        # print('self.edg[0] = ', self.edg[0], type(self.edg[0]))
        # print('self.edg[-1] = ', self.edg[-1], type(self.edg[-1]))
        # raise ValueError

        for i in range( len(self.edg)-2, -1, -1 ):
            print("#", str(i))
            for v in self.vex[i+1]:
                print('vex[i+1] = ', v)
            # print('self.edg[0] = ', self.edg[0], type(self.edg[0]))
            # print('self.edg[4] = ', self.edg[4], type(self.edg[4]))
            for e in self.edg[i+1]:
                print('edg[i+1] = ', e)
            vex_next = vex_of_edg( self.vex[i+1], self.edg[i+1])
            print('vex_next = ', v)
            for v in vex_next:
                print('vex_next = ', v)
            edg_new_list = []
            for edge in self.edg[i]:
                v0, a, v1 = edge
                if in_list( vex_next, v1 ):
                    edg_new_list.append(edge)
            # print('edg_new_list = ', edg_new_list)
            self.edg[i] = edg_new_list
        # for i, e_row in enumerate(self.edg):
        #     print("#e", i, " e = ", e_row)        
        # raise Exception

        for i, v_row in enumerate(self.vex):
            for v in v_row:
                print("#v", i, " v = ", v)
            for e in self.edg[i]:
                print("#e", i, " e = ", e)
