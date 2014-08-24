"""
Some support classes for easy array manipulation
"""

def _same_len(operation):
    """ Make sure all arrays are the same length """
    def decorator(func):
        def new_func(self, other):
            if len(self) != len(other):
                raise ValueError('array size mismatch for %s' % operation)
            return func(self, other)

        return new_func

    return decorator

class NumberArray(list):
   
    def __init__(self, stuff=None):
        if stuff is None:
            list.__init__(self)
        else:
            list.__init__(self, stuff)

    def __mul__(self, scalar):
        return NumberArray([i*scalar for i in self])

    __rmul__ = __mul__

    def __imul__(self, scalar):
        for i in xrange(len(self)): self[i] *= scalar
        return self
   
    @_same_len('addition')
    def __add__(self, other):
        return NumberArray([i+j for i,j in zip(self, other)])
   
    @_same_len('addition')
    def __iadd__(self, other):
        for i, val in enumerate(other): self[i] += val
        return self
   
    @_same_len('subtraction')
    def __sub__(self, other):
        return NumberArray([i-j for i,j in zip(self, other)])

    @_same_len('subtraction')
    def __isub__(self, other):
        for i, val in enumerate(other): self[i] -= val
        return self
