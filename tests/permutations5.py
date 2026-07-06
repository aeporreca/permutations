from sympy import *
from math import gcd, lcm
from frozendict import frozendict


def first(items):
    return items[0]


class Perm(AtomicExpr):

    is_commutative = True
    is_Number = True

    def __init__(self, mults={}):
        self.mults = frozendict(mults)

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
            if x.p == 0:
                return Perm()
            else:
                return Perm({1: x.p})
        raise TypeError(f'cannot convert {x!r} to {Perm}')

    def __repr__(self):
        if not self.mults:
            return '0'
        reprs = []
        for length, mult in sorted(self.mults.items(), key=first):
            if mult != 1 or length == 1:
                r = f'{mult!r}'
            else:
                r = ''
            if mult != 1 and length != 1:
                r += '*'
            if length != 1:
                r += f'C({length!r})'
            reprs.append(r)
        return ' + '.join(reprs)

    __str__ = __repr__

    def __hash__(self):
        return hash(self.mults)

    def _hashable_content(self):
        return tuple(sorted(self.mults.items()))

    def __add__(self, other):
        if isinstance(other, int | Integer):
            other = Perm.of(other)
        if not isinstance(other, Perm):
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


def C(n):
    if n == 0:
        return Perm()
    else:
        return Perm({n: 1})


X = Symbol('X')
Y = Symbol('Y')


# def tracefunc(frame, event, arg, indent=[0]):
#       if event == "call":
#           indent[0] += 2
#           print("-" * indent[0] + "> call function", frame.f_code.co_name)
#       elif event == "return":
#       #     print("<" + "-" * indent[0], "exit function", frame.f_code.co_name)
#           indent[0] -= 2
#       return tracefunc

# import sys
# sys.setprofile(tracefunc)


def patch_add(method):
    def wrapped(self, other):
        if isinstance(other, Perm):
            return other + self
        else:
            return method(self, other)
    return wrapped


def patch_mul(method):
    def wrapped(self, other):
        if isinstance(other, Perm):
            return other * self
        else:
            return method(self, other)
    return wrapped


Integer.__add__ = patch_add(Integer.__add__)
Integer.__mul__ = patch_mul(Integer.__mul__)
