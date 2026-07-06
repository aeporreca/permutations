from frozendict import frozendict
from itertools import product
from copy import deepcopy
from collections import defaultdict
from math import gcd, lcm


class Monomial:

    def __init__(self, exp={}):
        self.exp = frozendict(exp)

    @staticmethod
    def of(expr):
        if isinstance(expr, Monomial):
            return expr
        else:
            assert(expr == 1)
            return Monomial()

    def __hash__(self):
        return hash(self.exp)

    def __repr__(self):
        if not self.exp:
            return '1'
        reprs = []
        for var, exp in self.exp.items():
            if exp == 1:
                reprs.append(var)
            else:
                reprs.append(f'{var}**{exp}')
        return '*'.join(reprs)

    def __eq__(self, other):
        other = Monomial.of(other)
        return self.exp == other.exp

    def __mul__(self, other):
        other = Monomial.of(other)
        exp = {}
        for var in self.exp | other.exp:
            exp[var] = (self.exp.get(var, 0) +
                         other.exp.get(var, 0))
        return Monomial(exp)

    def vars(self):
        return sorted(self.exp.keys())

    def eval(self, value):
        res = 1
        for var in self.exp:
            res *= value[var] ** self.exp[var]
        return res


class Polynomial:

    def __init__(self, coeff=set()):
        self.coeff = frozendict(coeff)

    @staticmethod
    def of(x):
        if isinstance(x, Polynomial):
            return x
        else:
            return Polynomial({Monomial(): x})

    def __hash__(self):
        return hash(self.coeff)

    def __repr__(self):
        if not self.coeff:
            return '0'
        reprs = []
        for mono, coeff in self.coeff.items():
            if mono == 1:
                reprs.append(repr(coeff))
            elif coeff == 1:
                reprs.append(repr(mono))
            # elif len(coeff.cycles) > 1:
            #     reprs.append(f'({repr(coeff)})*{repr(mono)}')
            else:
                reprs.append(f'{repr(coeff)}*{repr(mono)}')
        return ' + '.join(reprs)

    def __add__(self, other):
        other = Polynomial.of(other)
        coeff = {}
        for mono in self.coeff | other.coeff:
            coeff[mono] = (self.coeff.get(mono, 0) +
                           other.coeff.get(mono, 0))
            if coeff[mono] == 0:
                del coeff[mono]
        return Polynomial(coeff)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-1) * other

    def __mul__(self, other):
        other = Polynomial.of(other)
        coeff = {}
        for mono1, mono2 in product(self.coeff, other.coeff):
            mono = mono1 * mono2
            coeff[mono] = (coeff.get(mono, 0) +
                           self.coeff.get(mono1, 1) *
                           other.coeff.get(mono2, 1))
            if coeff[mono] == 0:
                del coeff[mono]
        return Polynomial(coeff)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exp):
        res = 1
        for _ in range(exp):
            res *= self
        return res

    def vars(self):
        return sorted({var for mono in self.coeff for var in mono.vars()})

    def __call__(self, *values):
        if len(values) != len(self.vars()):
            raise ValueError('wrong number of arguments')
        values = {var: value for var, value in zip(self.vars(), values)}
        res = 0
        for mono in self.coeff:
            res += self.coeff[mono] * mono.eval(values)
        return res

    def __eq__(self, other):
        other = Polynomial.of(other)
        return self.coeff == other.coeff


class Permutation:

    def __init__(self, cycles={}):
        self.cycles = frozendict(cycles)

    @staticmethod
    def of(x):
        if isinstance(x, Permutation):
            return x
        else:
            return Permutation({1: x})

    def __repr__(self):
        if not self.cycles:
            return '0'
        elif list(self.cycles) == [1]:
            return repr(self.cycles[1])
        elif len(self.cycles) == 1:
            n = list(self.cycles)[0]
            if self.cycles[n] != 1:
                return f'{self.cycles[n]}*C({n})'
            else:
                return f'C({n})'
        else:
            return ' + '.join(repr(self.cycles[n] * Permutation({n: 1}))
                for n in reversed(sorted(self.cycles)))

    def __add__(self, other):
        other = Permutation.of(other)
        cycles = defaultdict(int, deepcopy(self.cycles))
        for n in other.cycles:
            cycles[n] += other.cycles[n]
            if cycles[n] == 0:
                del cycles[n]
        return Permutation(cycles)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-1) * other

    def __rsub__(self, other):
        return self - other

    def __mul__(self, other):
        other = Permutation.of(other)
        cycles = defaultdict(int)
        for m, n in product(self.cycles, other.cycles):
            cycles[lcm(m, n)] +=  self.cycles[m] * other.cycles[n] * gcd(m, n)
        return Permutation(cycles)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        other = Permutation.of(other)
        return self.cycles == other.cycles

    def __le__(self, other):
        other = Permutation.of(other)
        return all(self.cycles[length] <= other.cycles.get(length, 0)
                   for length in self.cycles)
        

def C(n):
    return Polynomial.of(Permutation({n: 1}))


# def extract_terms_with_cycle_of_length(P, length):
#     coeffs = {}
#     for mono in P.coeff:
#         coeff = Permutation.of(P.coeff[mono])
#         if length in coeff.cycles:
#             coeffs[mono] = coeff.cycles[length]
#     return Polynomial(coeffs)


def coefficient_of_vector(P, cycle):
    cycle = Permutation.of(cycle)
    # only works when the vector is a cycle
    assert (len(cycle.cycles) == 1 and
            list(cycle.cycles.values()) == [1])
    length = list(cycle.cycles.keys())[0]
    return Polynomial({mono: Permutation.of(coeff).cycles[length]
                       for mono, coeff in P.coeff.items()
                       if cycle <= Permutation.of(coeff)})
    

# def cycle(n):
#     return Permutation({n: 1})


def var(name):
    mono = Monomial({name: 1})
    return Polynomial({mono: 1})


X = var('X')
Y = var('Y')
Z = var('Z')

x1 = var('x1')
x2 = var('x2')
x3 = var('x3')
x6 = var('x6')
