import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

zero = np.zeros(2, dtype=int)

def vex_eq( v0, v1 ):
    print("Equal: ", v0, v1, np.array_equal(v0, v1))
    return np.array_equal(v0, v1)

def vex_of_edg( vex, edg ):
    edgs = []
    # print("vex = ", vex)
    for v in vex:
        for e in edg:
            v0, a, v1 = e
            # print("e = ", e)
            if vex_eq(v0,v):
                edgs.append(v)
                break
    return edgs

def array2str( in_arr, i ):
    name = str(i) + '_'
    for s in in_arr:
        name += str(s)
    return name

def arr2bistr( arr, dim ):
    bistr = ""
    for symbol in arr:
        for digit in symbol.value:
            bistr = bistr + str(int(digit))
    # print('bistr = ', bistr)
    return bistr

def in_list( vec_list, target ):
    for vec in vec_list:
        if np.array_equal(vec,target): return True
    return False

def read_mat( filename ):
    with open( filename, newline='') as f:
        reader = csv.reader(f)
        rows = []
        symbols = []
        for row in reader:
            rows = []
            for symbol in row:
                symbol_str = []
                for gf in symbol:
                    symbol_str.append(int(gf))
                print("symbol = ", symbol_str)
                rows.append(symbol_str)
            symbols.append(rows)
        print("symbols = ", symbols)
        rows_nparr = np.swapaxes(np.array(symbols),0,1)
        print("rows_nparr = ", rows_nparr.shape)
        size = int(int(rows_nparr.size)/len(symbols))
        print('size = ', size)
        print("p_mat = ", rows_nparr)
    return rows_nparr
