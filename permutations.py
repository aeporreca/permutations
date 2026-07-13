from frozendict import frozendict
from itertools import product, chain, combinations, repeat, count
from functools import total_ordering
from copy import deepcopy
from collections import defaultdict
from math import gcd, lcm, prod
import numpy as np
from scipy.optimize import milp, LinearConstraint
from sympy import Matrix


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

    def eval(self, values):
        res = 1
        for var in self.exp:
            res *= values.get(var, variable(var)) ** self.exp[var]
        return res

    def degree(self):
        return sum(self.exp.values())

    def is_linear(self):
        return self.degree() <= 1


class Polynomial:

    def __init__(self, coeff=set()):
        self.coeff = frozendict(coeff)

    @staticmethod
    def of(x):
        if isinstance(x, Polynomial):
            return x
        elif isinstance(x, Monomial):
            return Polynomial({x: 1})
        else:
            if x != 0:
                return Polynomial({Monomial(): x})
            else:
                return Polynomial()

    def as_monomial(self):
        error = ValueError(f'{self} is not a monomial')
        if len(self.coeff) != 1:
            raise error
        mono, coeff = list(self.coeff.items())[0]
        if coeff != 1:
            raise error
        return mono

    def __hash__(self):
        return hash(self.coeff)

    def __repr__(self):
        if not self.coeff:
            return '0'
        reprs = []
        for mono, coeff in self.coeff.items():
            coeff = Permutation.of(coeff)
            if mono == 1:
                reprs.append(repr(coeff))
            elif coeff == 1:
                reprs.append(repr(mono))
            # TODO: print negatives better
            elif len(coeff.cycles()) > 1:
                reprs.append(f'({repr(coeff)})*{repr(mono)}')
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
        return prod(repeat(self, exp))

    def __call__(self, *values, **kwargs):
        error = TypeError('wrong arguments')
        if values and kwargs:
            raise error
        if values:
            if len(values) != len(self.vars()):
                raise error
            values = {var: value for var, value in zip(self.vars(), values)}
        else:
            values = kwargs
        res = 0
        for mono in self.coeff:
            res += self.coeff[mono] * mono.eval(values)
        return res

    def __eq__(self, other):
        other = Polynomial.of(other)
        return self.coeff == other.coeff

    def vars(self):
        return sorted({var for mono in self.coeff
                       for var in mono.vars()})

    def is_linear(self):
        return all(mono.is_linear() for mono in self.coeff)

    def is_univariate(self):
        return len(self.vars()) <= 1

    def coefficients(self):
        return set(self.coeff.values())

    def terms(self):
        return self.coeff

    def coeff_of(self, mono):
        if isinstance(mono, Polynomial):
            mono = mono.as_monomial()
        elif isinstance(mono, int):
            mono = Monomial.of(mono)
        return self.coeff.get(mono, 0)

    def homomorphism(self, hom):
        return Polynomial({mono: hom(coeff)
                           for mono, coeff in self.coeff.items()})

    def degree(self):
        return max(mono.degree() for mono in self.coeff)

    def minimal_degree(self):
        return min(mono.degree() for mono in self.coeff)

    @staticmethod
    def terms_with_cycle(P, cycle):
        P = Polynomial.of(P)
        cycle = Permutation.of(cycle)
        if len(cycle.cycles()) != 1:
            raise ValueError(f'{cycle} is not a cycle')
        length = cycle.lengths()[0]
        res = 0
        for mono, coeff in P.terms().items():
            coeff = Permutation.of(coeff)
            if length in coeff.lengths():
                res += coeff.multiplicity(length) * Polynomial.of(mono)
        return res


def partitions(n, m=1):
    if n == 0:
        yield ()
    else:
        for k in range(m, n + 1):
            for p in partitions(n - k, k):
                yield (k,) + p


@total_ordering
class Permutation:

    def __init__(self, cycles={}):
        self._cycles = frozendict(cycles)

    @staticmethod
    def generate(size):
        for p in partitions(size):
            yield sum(C(n) for n in p)

    @staticmethod
    def of(x):
        if isinstance(x, Permutation):
            return x
        elif isinstance(x, int):
            if x != 0:
                return Permutation({1: x})
            else:
                return Permutation()
        else:
            return NotImplemented

    def __repr__(self):
        if not self._cycles:
            return '0'
        elif list(self._cycles) == [1]:
            return repr(self._cycles[1])
        elif len(self._cycles) == 1:
            n = list(self._cycles)[0]
            if self._cycles[n] != 1:
                return f'{self._cycles[n]}*C({n})'
            else:
                return f'C({n})'
        else:
            return ' + '.join(repr(self._cycles[n] * Permutation({n: 1}))
                for n in reversed(sorted(self._cycles)))

    def __hash__(self):
        return hash(self._cycles)

    def __len__(self):
        return sum(mult * length
                   for mult, length in self._cycles.items())

    def __add__(self, other):
        other = Permutation.of(other)
        if other is NotImplemented:
            return NotImplemented
        cycles = defaultdict(int, deepcopy(self._cycles))
        for n in other._cycles:
            cycles[n] += other._cycles[n]
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
        if other is NotImplemented:
            return NotImplemented
        cycles = defaultdict(int)
        for m, n in product(self._cycles, other._cycles):
            l = lcm(m, n)
            cycles[l] += (self._cycles[m] *
                          other._cycles[n] * gcd(m, n))
            if cycles[l] == 0:
                del cycles[l]
        return Permutation(cycles)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exp):
        return prod(repeat(self, exp))

    def __eq__(self, other):
        other = Permutation.of(other)
        return self._cycles == other._cycles

    def __le__(self, other):
        other = Permutation.of(other)
        cycles = sorted(self._cycles | other._cycles)
        mults1 = [self._cycles.get(length, 0) for length in cycles]
        mults2 = [other._cycles.get(length, 0) for length in cycles]
        return mults1 >= mults2

    def cycles(self):
        return frozenset(C(n) for n in self._cycles)

    def lengths(self):
        return sorted(self._cycles.keys())

    def multiplicity(self, length):
        return self._cycles.get(length, 0)

    def sqrt(self, k=2):
        res = Permutation.of(0)
        while len(res**k) < len(self):
            length = (self - res**k).lengths()[0]
            res += C(length)
        if res**k != self:
            return None
        return res

    def minimal_equation(self):
        for deg in count(1):
            E = NaturalExtension(self.cycles())
            B = sorted(E.basis())
            dim = len(B)
            nvars = deg + 1
            A = Matrix.zeros(dim, nvars, dtype=int)
            for j in range(nvars):
                powA = Permutation.of(self**j)
                for i in range(dim):
                    A[i, j] = powA.multiplicity(len(B[i]))
            sol = A.nullspace()
            assert len(sol) <= 1
            if not sol:
                continue
            a = list(map(int, sol[0]))
            P = sum(a[i]*X**i for i in range(deg + 1) if a[i] > 0)
            P = Polynomial.of(P)
            Q = sum(-a[i]*X**i for i in range(deg + 1) if a[i] < 0)
            Q = Polynomial.of(Q)
            assert P(**{'X': self}) == Q(**{'X': self})
            return Equation(P, Q)


def C(n):
    if n == 0:
        return Permutation()
    else:
        return Permutation({n: 1})


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def prime_power_factors(n):
    divisors = [1]
    for d in range(2, n + 1):
        e = 1
        while n % (e * d) == 0:
            e = e * d
        if e != 1:
            divisors.append(e)
        n = n // e
        if n == 1:
            break
    return divisors


def divisors_of_cycles(perms):
    cycles = set(chain.from_iterable(
        Permutation.of(perm).cycles() for perm in perms))
    divs = cycles.copy()
    for cycle in cycles:
        divs |= {C(n) for n in prime_power_factors(cycle.lengths()[0])}
    return divs


class NaturalExtension:

    def __init__(self, perms):
        cycles = chain.from_iterable(perm.cycles() for perm in perms)
        self.generators = frozenset(cycles)

    def __repr__(self):
        return f'NaturalExtension({set(self.generators)})'

    def basis(self):
        return set(chain.from_iterable(
            Permutation.of(prod(X)).cycles()
            for X in powerset(self.generators)))

    def degree(self):
        return len(self.basis())


def variable(name):
    mono = Monomial({name: 1})
    return Polynomial({mono: 1})


class Equation:

    def __init__(self, P, Q):
        self.P = Polynomial.of(P)
        self.Q = Polynomial.of(Q)
        D = divisors_of_cycles(self.P.coefficients() |
                               self.Q.coefficients())
        E = NaturalExtension(D)
        self.B = sorted(E.basis())

    def __repr__(self):
        return f'Equation({self.P}, {self.Q})'

    def is_linear(self):
        return self.P.is_linear() and self.Q.is_linear()

    def is_univariate(self):
        return len(set(self.P.vars() + self.Q.vars())) <= 1

    def vars(self):
        return sorted(set(self.P.vars() + self.Q.vars()))

    def natural_variables(self):
        vars = self.vars()
        return {(var, cycle): variable(f'V[{var},{cycle}]')
                for var in vars for cycle in self.B}

    def degree(self):
        return max(self.P.degree(), self.Q.degree())

    def as_natural_equations(self):
        V = self.natural_variables()
        W = {var: sum(cycle * V[var,cycle] for cycle in self.B)
             for var in self.vars()}
        P1 = self.P(*(W[var] for var in self.P.vars()))
        Q1 = self.Q(*(W[var] for var in self.Q.vars()))
        equations = []
        for cycle in self.B:
            P2 = Polynomial.terms_with_cycle(P1, cycle)
            Q2 = Polynomial.terms_with_cycle(Q1, cycle)
            equations.append((Polynomial.of(P2), Polynomial.of(Q2)))
        return equations

    def as_matrix_equation(self):
        dim = len(self.B)
        nvars = len(self.vars())
        A = np.zeros((dim, dim * nvars), dtype=int)
        a = np.zeros(dim, dtype=int)
        B = np.zeros((dim, dim * nvars), dtype=int)
        b = np.zeros(dim, dtype=int)
        V = self.natural_variables()
        for i, (P, Q) in enumerate(self.as_natural_equations()):
            for j, var in enumerate(self.natural_variables()):
                A[i][j] = P.coeff_of(V[var])
                B[i][j] = Q.coeff_of(V[var])
            a[i] = P.coeff_of(1)
            b[i] = Q.coeff_of(1)
        return A - B, b - a

    def solve_linear(self):
        A, b = self.as_matrix_equation()
        dim, nvars = A.shape
        constraints = LinearConstraint(A, b, b)
        objective = np.zeros(nvars, dtype=int) # trivial
        integrality = np.ones(nvars, dtype=int)
        res = milp(c=objective,
                   constraints=constraints,
                   integrality=integrality)
        if res.x is None:
            return None
        vars = self.vars()
        x = np.array(res.x, dtype=int)
        sol = {vars[i]: x[i*dim:(i+1)*dim] @ self.B
               for i in range(len(self.vars()))}
        assert self.P(**sol) == self.Q(**sol)
        yield sol

    def solutions_univariate(self):
        if not self.is_univariate():
            raise ValueError(f'{self} is not univariate')
        P = self.P.homomorphism(card)
        Q = self.Q.homomorphism(card)
        var = P.vars()[0]
        for size in solutions_natural_univariate_equation(P, Q):
            for perm in Permutation.generate(size):
                if self.P(**{var: perm}) == self.Q(**{var: perm}):
                    yield perm

    # def solve(self):
    #     if self.is_univariate():
    #         return self.solve_univariate()
    #     elif self.is_linear():
    #         return self.solve_linear()
    #     raise NotImplementedError(
    #         f'unable to solve {self}')


def card(p):
    if isinstance(p, int):
        return p
    elif isinstance(p, Permutation):
        return len(p)


def divisors(n):
    for d in range(1, n + 1):
        if n % d == 0:
            yield d
    

def solutions_natural_univariate_equation(P, Q):
    R = P - Q
    if R == 0:
        yield from count(0)
    a = R.coeff_of(1)
    if a == 0:
        yield 0
        var = R.vars()[0]
        tmp = {}
        k = R.minimal_degree()
        for mono, coeff in R.coeff.items():
            if mono.exp[var] > k:
                tmp[Monomial({var: mono.exp[var] - k})] = coeff
            else:
                tmp[Monomial()] = coeff
        S = Polynomial(tmp)
        yield from solutions_natural_univariate_equation(S, 0)
    else:
        for d in divisors(abs(a)):
            if R(d) == 0:
                yield d


def solutions_natural_univariate_system(equations):
    if not equations:
        yield {}
    else:
        P, Q = equations[0]
        vars = sorted(set(P.vars() + Q.vars()))
        if not vars:
            raise ValueError(f'underdetermined system {equations}')
        var = vars[0]
        for x in solutions_natural_univariate_equation(P, Q):
            equations2 = [(Polynomial.of(P1(**{var: x})),
                           Polynomial.of(Q1(**{var: x})))
                        for P1, Q1 in equations[1:]]
            for sol in solutions_natural_univariate_system(equations2):
                yield sol | {var: x}


X = variable('X')
Y = variable('Y')
Z = variable('Z')
