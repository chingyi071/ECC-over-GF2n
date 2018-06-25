import numpy as np
from util import *

class GFn:
    def __init__( self, value, nbit ):
        self.nbit = nbit
        # print('value = ', value)
        if len(np.shape(value)) is not 1:
            raise ValueError()
        self.value = np.array(value).flatten().astype(int)

    def __add__( self, a ):
        if self.value.shape[0] is not a.value.shape[0]:
            err_msg = "Mismatched size in GFn addition"
            err_msg += ", augend is ", self.value.shape
            err_msg += ", addend is ", a.value.shape
            raise ValueError()
        result_value = fit_gfn( a.value + self.value, self.nbit)
        result = GFn( result_value, self.nbit)
        return result

    def __mul__( self, a ):

        if len(self.value.shape) > 1:
            raise ValueError

        # If multiplicend is an array, separte it into elements
        if type(a) is type(np.array([])):
            flat_input = a.flatten()

            # flat_output is an 1-D array storing each element result
            flat_output = np.empty_like( a, dtype=object )
            for i in range(0, flat_output.size):
                flat_output[i] = self*flat_input[i]

            # Reshape the result as the shape of the input
            return np.reshape( flat_output, a.shape )

        if type(a) is not type(GFn([0],1)):
            raise ValueError( str(a) + "is neither GFn or np.array")

        # Shift multiplicend and add it to psum
        product = np.zeros_like(a.value, dtype=int).astype(int)
        for i in range( 0, self.nbit ):

            # Partial sum = [0]*i + self*a
            psum = np.append( np.zeros(i), self.value[i]*a.value)
            psum = fit_gfn( psum, self.nbit )
            product = product + psum

        if product.shape[0] is not a.nbit:
            raise Exception()

        return GFn( product, a.nbit )

    def __eq__( self, a ):
        if not self.nbit == a.nbit:
            raise ValueError("GFn compare two input with different bit length")
            return False
        return np.array_equal(self.value, a.value)

    def __repr__( self ):
        return 'GFn(' + str(self) + ')'

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
