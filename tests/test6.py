from sympy import *
from math import gcd, lcm
from functools import total_ordering

@total_ordering
class C(Symbol):

    def __new__(cls, length):
        c = super().__new__(cls, f'C({length})')
        c.length = length
        return c

    # @staticmethod
    # def of(x):
    #     if isinstance(x, C):
    #         return x
    #     elif isinstance(x, int):
    #         return C(x)
    #     elif isinstance(x, Integer):
    #         return C(x.p)
    #     raise TypeError(f'cannot convert {x!r} to {C}')

    # def __hash__(self):
    #     return hash(self.length)

    # def _hashable_content(self):
    #     return (self.length, )

    # def __add__(self, other):
    #     return Add(self, other)

    # def __add__(self, other):
    #     print(f'__add__({self}, {other})')
    #     if isinstance(other, int | Integer):
    #         other = C.of(other)
    #     if not isinstance(other, C):
    #         # return other + self
    #         return NotImplemented
    #     mults = {}
    #     for length in self.mults | other.mults:
    #         mults[length] = (self.mults.get(length, 0) +
    #                          other.mults.get(length, 0))
    #     return C(mults)

    # def __radd__(self, other):
    #     return self + other
    
    def __hash__(self):
        return hash(self.length)

    def __mul__(self, other):
        if isinstance(other, C):
            mult = gcd(self.length, other.length)
            length = lcm(self.length, other.length)
            return mult * C(length)
        # elif isinstance(other, Mul):
        else:
            return super().__mul__(other)

    def __rmul__(self, other):
        return self * other

    def __lt__(self, other):
        if not isinstance(other, C):
            return False
        return self.length < other.length

    def __eq__(self, other):
        if not isinstance(other, C):
            return False
        return self.length == other.length


def demulify(mul):
    if not isinstance(mul, Mul):
        return mul
    if len(mul.args) == 2:
        if isinstance(mul.args[0], Integer):
            return mul
    # print(mul.args[0] * mul.args[1])
    if len(mul.args) == 2:
        return mul.args[0] * mul.args[1]
    else:
        return Mul(mul.args[0] * mul.args[1], demulify(*mul.args[2:]))
