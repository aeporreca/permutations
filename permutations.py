from math import gcd, lcm
from collections import defaultdict
from itertools import product, count
from functools import total_ordering
from copy import copy
from sys import stderr


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
        r = C(1)
        for i in range(n):
            r = r * self
        return r

    def __eq__(self, other):
        if isinstance(other, int):
            other = other * C(1)
        return self.cycles == other.cycles

    def __le__(self, other):
        for n in self.cycles:
            if not self.cycles[n] <= other.cycles[n]:
                return False
        return True

    def __len__(self):
        return sum(n * self.cycles[n] for n in self.cycles)

    def coeff(self, n):
        return self.cycles[n]

    def lengths(self):
        return sorted(self.cycles.keys())


class Poly:

    def __init__(self, *coeff):
        self.coeff = coeff

    def __call__(self, val):
        res = 0
        n = len(self.coeff)
        for i in range(n):
            res = res * val + self.coeff[i]
        return res

    def __repr__(self):
        return f'Poly{self.coeff}'


def minimal_poly(elem, degree):
    for n in count(0):
        print(f'{n = }', file=stderr)
        for coeff in product(range(-n + 1, n), repeat=degree):
            P = Poly(1, *coeff)
            if P(elem) == 0:
                return P


def linear_combination(elem, degree):
    n = list(sorted(elem.cycles))[0]
    bound = (n * elem.cycles[n])**(degree - 1)
    for n in range(bound + 1):
        print(f'{n = }', file=stderr)
        for coeff in product(range(-n + 1, 1), repeat=degree):
            P = Poly(1, *coeff)
            if P(elem) == 0:
                return P


# Minimal polynomial of a given degree

import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint

def minimal_poly(elem, deg):
    rows = (elem**deg).lengths()[-1] + 1
    A = np.zeros((rows + 1, deg + 1))
    for i in range(rows):
        for j in range(deg + 1):
            A[i][j] = (elem**j).coeff(i)
    A[rows][-1] = 1
    L = np.zeros(rows + 1)
    U = np.zeros(rows + 1)
    L[-1] = 1
    U[-1] = np.inf
    constraints = LinearConstraint(A, L, U)
    c = np.zeros(deg + 1)
    integrality = np.ones(deg + 1)
    res = milp(c=c, bounds=Bounds(),
               constraints=constraints,
               integrality=integrality)
    return Poly(*(int(n) for n in reversed(res.x)))
