from math import gcd, lcm
from collections import defaultdict
from itertools import product
from functools import total_ordering
from copy import copy


@total_ordering
class C:

    def __init__(self, n=0):
        self.cycles = defaultdict(int)
        if n > 0:
            self.cycles[n] = 1
        self.normalize()

    def __repr__(self):
        if not self.cycles:
            return '0'
        elif list(self.cycles) == [1]:
            return repr(self.cycles[1])
        elif len(self.cycles) == 1:
            n = list(self.cycles)[0]
            if self.cycles[n] != 1:
                return f'{self.cycles[n]} * C({n})'
            else:
                return f'C({n})'
        else:
            return ' + '.join(repr(self.cycles[n] * C(n))
                for n in reversed(sorted(self.cycles)))

    def normalize(self):
        for n in list(self.cycles):
            if self.cycles[n] == 0:
                del self.cycles[n]

    def __add__(self, other):
        if isinstance(other, int):
            other = other * C(1)
        perm = copy(self)
        for n in other.cycles:
            perm.cycles[n] += other.cycles[n]
        perm.normalize()
        return perm

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + (-1) * other

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other, int):
            tmp = C(1)
            tmp.cycles[1] = other
            tmp.normalize()
            other = tmp
        perm = C()
        for m, n in product(self.cycles, other.cycles):
            perm.cycles[lcm(m, n)] += (
                self.cycles[m] * other.cycles[n] * gcd(m, n)
            )
        perm.normalize()
        return perm

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        r = 1
        for i in range(n):
            r = r * self
        return r

    def __eq__(self, other):
        return self.cycles == other.cycles

    def __le__(self, other):
        for n in self.cycles:
            if not self.cycles[n] <= other.cycles[n]:
                return False
        return True
