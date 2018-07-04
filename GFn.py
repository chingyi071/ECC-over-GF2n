import numpy as np
from util import *

class GFn:
    def __init__( self, value, nbit, dump=False ):
        self.nbit = nbit
        # print('value = ', value)
        if np.isscalar(value):
            bin_list = [int(i) for i in bin(value)[2:]]
            value_array = np.flip( [0]*(nbit-len(bin_list)) + bin_list, axis=0 )

        elif len(np.shape(value)) is 1:
            if np.shape(value)[0] > nbit:
                value_array = fit_gfn( value, nbit )
            elif np.shape(value)[0] < nbit:
                value_array = np.append( value, np.zeros(nbit-np.shape(value)[0]) )
            else:
                value_array = np.array(value)

        if value_array.shape[0] is not nbit:
            print("shape = ", value_array.shape)
            print("ideal = ", (nbit,))
            raise Exception
        self.value = np.array(value_array).astype(int)
        self.dump = dump

    def __radd__( self, a ):
        if self.dump:
            print("radd(", str(self), ",", a, type(a), ")")
        if int(a) == 0: return GFn( self.value, self.nbit)
        print("radd: Add with non-zero")
        raise Exception()

    def __add__( self, a ):
        if self.dump:
            print("add(", str(self), self.value.shape, ",", a, type(a), ")")
        if type(a) is not type(GFn(np.array([0]),1)):
            if int(a) == 0: return GFn(self.value, self.nbit)        
        if self.value.shape[0] is not a.value.shape[0]:
            err_msg = "Mismatched size in GFn addition"
            err_msg += ", augend is ", self.value.shape
            err_msg += ", addend is ", a.value.shape
            raise ValueError()
        result_value = fit_gfn( a.value + self.value, self.nbit)
        result = GFn( result_value, self.nbit)
        if self.dump:
            print("add return ", str(result))
        return result

    def __sub__( self, a ):
        return self+a

    def __rmul__( self, a ):
        if self.dump:
            print("rmul(", str(self), ",", a, type(a), ")")
        result_value = np.remainder( a*self.value, 2 )
        result = GFn( result_value, self.nbit )
        if self.dump:
            print("rmul return ", str(result), type(result))
        return result 

    def __mul__( self, a ):
        if self.dump:
            print("mul(", str(self), self.value.shape, ",", a, type(a), ")")
        if len(self.value.shape) > 1:
            raise ValueError

        # If multiplicend is an array, separte it into elements
        if type(a) is type(np.array([])):
            return a.__rmul__(self)

        if np.isscalar(a):
            value = a*self.value
            return GFn( fit_gfn(value, self.nbit), self.nbit)

        if type(a) is not type(GFn(0,1)):
            raise ValueError( str(a) + "is neither GFn or np.array")

        # Shift multiplicend and add it to psum
        product = np.zeros_like(self.value, dtype=int).astype(int)
        for i in range( 0, self.nbit ):

            # Partial sum = [0]*i + self*a
            psum = np.append( np.zeros(i), self.value[i]*a.value)
            psum = fit_gfn( psum, self.nbit )
            if psum.size != self.nbit:
                print("Size mismatch1")
                raise Exception
            product = product + psum
            product = fit_gfn( product, self.nbit )
            if product.size != self.nbit:
                print("Size mismatch2")
                raise Exception
        product = fit_gfn( product, self.nbit )

        if product.shape[0] is not a.nbit:
            print("mul(", str(self), self.value.shape, ",", a, type(a), ")")
            print("self.nbit = ", self.nbit)
            print("a.nbit = ", a.nbit)
            print("Size mismatch3")
            raise Exception()

        result = GFn( product, a.nbit )
        if self.dump:
            print("Return GFn(", str(result), ")")

        return GFn( product, a.nbit )

    def __eq__( self, a ):
        if type(a) is not type(GFn(0,1)):
            return int(self) == a
        if not self.nbit == a.nbit:
            raise ValueError("GFn compare two input with different bit length")
            return False
        return np.array_equal(self.value, a.value)

    def __repr__( self ):
        return 'GFn(' + str(self) + ')'

    def __float__(self):
        return float(int(self))

    def __int__(self):
        total = 0
        for i, digit in enumerate(reversed(self.value)):
            total = total*2 + int(digit)
        return total

    def __str__(self):
        out_str = ""
        for x in self.value:
            out_str = str(x) + out_str
        return out_str

    def iszero(self):
        if np.count_nonzero(self.value) > 0: return False
        return True

    def toGF2(self,n):
        if n < self.nbit:
            return GFn( self.value[:n], n )
        else:
            return GFn( np.append(self.value, np.zeros(n-self.nbit)), n )

    def power(self,n):
        e = GFn(1,self.nbit)
        for i in range(0,n):
            e = e*self
        return e

    def __exp__(self,n):
        return self.power(n)

    def is_root(self, g):
        y = np.polyval(g,self)
        return y.iszero()

    def inverse(self):
        for i in range(0,2**self.nbit):
            if int( GFn(i,self.nbit) * self )==1:
                return GFn(i,self.nbit)

def intlist_to_gfpolylist( int_list, m ):
    return [GFn(g,m) for g in int_list]

def symbol_all( nbit ):
    ret_list = []
    for i in range(0,2**nbit):
        ret_list.append(GFn(i,nbit))
    return ret_list

def find_characteristic( a ):
    i = 1
    product = a
    while 1:
        if int(product) == 1: return i+1
        product = product*a
        i = i+1

def gfn_array_modulo( dividend, modular_poly ):
    logq = dividend[0].nbit
    zero_logq = GFn(0,logq)
    one_logq  = GFn(1,logq)

    modular = np.array(modular_poly)
    while 1:

        # Remainder is 0, return with padding or slicing
        if np.argwhere(dividend!=zero_logq).size == 0:
            if len(dividend) < len(modular):
                return np.append( [zero_logq] * (len(modular)-len(dividend)-1), dividend )
            else:
                return dividend[-len(modular):]

        msb = np.min(np.argwhere(dividend!=zero_logq),axis=0)[0]

        # Degree of modular is less than dividend
        if msb > len(dividend) - len(modular):
            if len(dividend) < len(modular):
                return np.append( [zero_logq] * (len(modular)-len(dividend)-1), dividend )
            else:
                return dividend[-len(modular):]

        # Do padding to align the MSB of modular to the dividend
        remainder = np.append( [zero_logq] * msb, modular )
        remainder = np.append( remainder, [zero_logq] * (len(dividend) - msb - len(modular)))
        remainder = remainder * dividend[msb]

        # Obtain the result from polynomial addition
        result = np.empty(len(dividend), dtype=object)
        for i, x in enumerate(result):
            result[i] = dividend[i]+remainder[i]
        dividend = result

def gen_zero_one_alpha_overGFq( q ):
    import math
    logq = int(math.log2(q))
    if q<4:
        return GFn(0,logq), GFn(1,logq), None
    else:
        return GFn(0,logq), GFn(1,logq), GFn(2,logq)

def gf_map( a, b, verbose=0 ):

    if a>b:
        s = int((2**a-1)/(2**b-1))
        if b==1:
            src = [GFn(1,a)]
            trg = [GFn(1,1)]
        else:
            alpha = GFn(2,a)
            src = [alpha.power(i) for i in range(0,2**a-1,s)]
            trg = [GFn(2,b).power(i) for i in range(0,2**b-1) ]
    elif a<b:
        s = int((2**b-1)/(2**a-1))
        if a==1:
            src = [GFn(1,1)]
            trg = [GFn(1,b)]
        else:
            alpha = GFn(2,a)
            src = [alpha.power(i) for i in range(0,2**a-1)]
            trg = [GFn(2,b).power(i) for i in range(0,2**b-1,s) ]
    else: raise ValueError

    if verbose:
        print("alpha^"+str(s), "in GF( 2^"+str(a),") = beta in GF( 2^"+str(b), ")")
    table = [(GFn(0,a), GFn(0,b) )]
    for i,xy in enumerate(zip(src,trg)):
        x,y = xy
        if verbose:
            print("#", i, "(alpha^"+str(s)+")^"+str(i), " = ", x, "on GF( 2^"+str(a),") = beta^"+str(i), "in GF( 2^"+str(b), ")")
        table.append((x,y))
    if verbose:
        print("--")
    return table

def weight( f ):
    int_list = [int(x) for x in f]
    non_zero_list = [x>0 for x in int_list]
    return sum(non_zero_list)

