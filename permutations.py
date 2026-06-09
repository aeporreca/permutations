from math import gcd, lcm
from collections import defaultdict
from itertools import product
from copy import copy


class C:

    def __init__(self, n=0):
        self.cycles = defaultdict(int)
        if n > 0:
            self.cycles[n] = 1
        self.normalize()

    def __repr__(self):
        if not self.cycles:
            return '0'
        elif list(self.cycles.keys()) == [1]:
            return repr(self.cycles[1])
        elif len(self.cycles.keys()) == 1:
            k = list(self.cycles.keys())[0]
            if self.cycles[k] != 1:
                return f'{self.cycles[k]} * C({k})'
            else:
                return f'C({k})'
        else:
            return ' + '.join(repr(self.cycles[k] * C(k))
                for k in reversed(sorted(self.cycles.keys())))

    def normalize(self):
        for k in list(self.cycles.keys()):
            if self.cycles[k] == 0:
                del self.cycles[k]

    def __add__(self, other):
        if isinstance(other, int):
            other = other * C(1)
        p = copy(self)
        for k in other.cycles:
            p.cycles[k] += other.cycles[k]
        p.normalize()
        return p

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if isinstance(other, int):
            tmp = C(1)
            tmp.cycles[1] = other
            other = tmp
            other.normalize()
        p = C()
        for k, h in product(self.cycles.keys(), other.cycles.keys()):
            p.cycles[lcm(k, h)] += self.cycles[k] * other.cycles[h] * gcd(k, h)
        p.normalize()
        return p

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        r = 1
        for i in range(n):
            r = r * self
        return r
