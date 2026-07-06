from sympy import *
from functools import total_ordering
# from frozendict import frozendict

# @total_ordering
class Perm(Expr):

    is_commutative = True
    is_number = False

    def __new__(cls, mults={}):
        return super().__new__(cls, frozendict(mults))

    @property
    def mults(self):
        return self.args[0]

    def __str__(self):
        return 'P'

    @staticmethod
    def of(x):
        if isinstance(x, Perm):
            return x
        elif isinstance(x, int):
            if x == 0:
                return Perm()
            else:
                return Perm({1: x})
        elif isinstance(x, Integer):
            return Perm({1: x.p})
        raise TypeError(f'cannot convert {x!r} to {cls}')

    def __add__(self, other):
        if isinstance(other, int | Integer):
            other = Perm.of(other)
        if not isinstance(other, Perm):
            # return other + self
            return NotImplemented
        mults = {}
        for length in self.mults | other.mults:
            mults[length] = (self.mults.get(length, 0) +
                             other.mults.get(length, 0))
        return Perm(mults)

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        if isinstance(other, int | Integer):
            other = Perm.of(other)
        if not isinstance(other, Perm):
            return NotImplemented
        mults = {}
        for length1, mult1 in self.mults.items():
            for length2, mult2 in other.mults.items():
                length = lcm(length1, length2)
                mult = gcd(length1, length2)
                mults[length] = (mults.get(length, 0) +
                                 mult * mult1 * mult2)
        return Perm(mults)

    def __rmul__(self, other):
        return self * other

    def __lt__(self, other):
        if isinstance(other, int | Integer):
            other = Perm.of(other)
        if not isinstance(other, Perm):
            return True
        return self._hashable_content() < other._hashable_content()

    def __eq__(self, other):
        return self.mults == other.mults

    def _hashable_content(self):
        return tuple(sorted(self.mults.items()))

    def __hash__(self):
        return hash(self._hashable_content())
