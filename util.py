import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

zero = np.zeros(2, dtype=int)

def vex_in_edg( vex, edg ):
    edgs = []
    for v in vex:
        for e in edg:
            v0, a, v1 = e
            if np.array_equal(v0,v):
                edgs.append(v)
                break
    return edgs

def array2str( in_arr, i ):
    return str(i) + '_' + str(in_arr[0]) + str(in_arr[1])

def int2bistr( i, dim ):
    bistr = ""
    input_int = i
    for nbit in range(0,dim):
        bistr = str(i%2) + bistr
        i = int(i/2)
    return bistr

def arr2bistr( arr, dim ):
    bistr = ""
    for digit in arr:
        bistr = bistr + str(digit)
    return bistr

def in_list( vec_list, target ):
    for vec in vec_list:
        if np.array_equal(vec,target): return True
    return False

def read_mat( filename ):
    with open( filename, newline='') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append([int(x) for x in row])
        p_mat = np.array(rows)
    return p_mat
