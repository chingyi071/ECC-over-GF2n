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
            value_array = value

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
            flat_input = a.flatten()

            # flat_output is an 1-D array storing each element result
            flat_output = np.empty_like( a, dtype=object )
            if flat_output.size == 1:
                flat_output = self*flat_input[0]
            else:
                for i in range(0, flat_output.size):
                    flat_output[i] = self*flat_input[i]

            # Reshape the result as the shape of the input
            return np.reshape( flat_output, a.shape )

        if np.isscalar(a):
            value = a*self.value
            return GFn( fit_gfn(value, self.nbit), self.nbit)

        if type(a) is not type(GFn([0],1)):
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
            ggg
            raise Exception()

        result = GFn( product, a.nbit )
        if self.dump:
            print("Return GFn(", str(result), ")")

        return GFn( product, a.nbit )

    def __eq__( self, a ):
        if type(a) is not type(GFn([0],1)):
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
        return GFn( self.value[:n], n )

def intlist_to_gfpoly( int_list, m ):
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
