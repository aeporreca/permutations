from math import gcd, lcm
from collections import defaultdict
from frozendict import frozendict
from itertools import product, chain
from functools import total_ordering
from copy import deepcopy


def ensure_perm(other):
    if isinstance(other, int):
        tmp = Perm(1)
        tmp.cycles[1] = other
        tmp.normalize()
        return tmp
    else:
        return other


@total_ordering
class Perm:

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
            return ' + '.join(repr(self.cycles[n] * Perm(n))
                for n in reversed(sorted(self.cycles)))

    def normalize(self):
        for n in list(self.cycles):
            if self.cycles[n] == 0:
                del self.cycles[n]

    def __add__(self, other):
        other = ensure_perm(other)
        perm = deepcopy(self)
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
        other = ensure_perm(other)
        perm = Perm()
        for m, n in product(self.cycles, other.cycles):
            perm.cycles[lcm(m, n)] += (
                self.cycles[m] * other.cycles[n] * gcd(m, n)
            )
        perm.normalize()
        return perm

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        r = Perm(1)
        for i in range(n):
            r = r * self
        return r

    def __eq__(self, other):
        if isinstance(other, int):
            other = other * Perm(1)
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

    def __hash__(self):
        return hash(frozendict(self.cycles))


def C(n):
    return ensure_poly(Perm(n))


def term_repr(coeff, monomial):
    if coeff == 1 and not monomial:
        return '1'
    s = ''
    if coeff != 1:
        s = repr(coeff)
        if monomial:
            s += ' * '
    vars = []
    for var in monomial:
        t = var
        if monomial[var] != 1:
            t += '**' + repr(monomial[var])
        vars.append(t)
    return s + ' * '.join(vars)


def ensure_poly(other):
    if isinstance(other, int | Perm):
        return Poly(defaultdict(Perm, {frozendict(): ensure_perm(other)}))
    else:
        return other


class Poly:

    def __init__(self, monomials=None):
        if monomials is None:
            monomials = defaultdict(Perm)
        for monomial in list(monomials.keys()):
            if monomials[monomial] == 0:
                del monomials[monomial]
        self.monomials = monomials

    def __repr__(self):
        if not self.monomials:
            return '0'
        terms = (term_repr(self.monomials[monomial], monomial)
                 for monomial in self.monomials)
        return ' + '.join(terms)

    def __add__(self, other):
        other = ensure_poly(other)
        res = deepcopy(self)
        for monomial in other.monomials:
            res.monomials[monomial] += other.monomials[monomial]
            if res.monomials[monomial] == 0:
                del res.monomials[monomial]
        return res

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-1) * other

    def __rsub__(self, other):
        return self - other

    def __mul__(self, other):
        other = ensure_poly(other)
        monomials = defaultdict(Perm)
        for mono1, mono2 in product(self.monomials, other.monomials):
            monomial = defaultdict(int, mono1)
            for var in mono2:
                monomial[var] += mono2[var]
            monomials[frozendict(monomial)] = self.monomials[mono1] * other.monomials[mono2]
        return Poly(monomials)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exp):
        res = ensure_poly(1)
        for i in range(exp):
            res *= self
        return res

    def vars(self):
        return sorted({var for monomial in self.monomials
                       for var in monomial})

    def coeffs(self):
        return set(self.monomials.values())

    def lengths(self):
        return list(chain.from_iterable(coeff.lengths() for coeff in self.coeffs()))
        
    def __call__(self, *values):
        vars = self.vars()
        if len(vars) != len(values):
            raise ValueError('wrong number of variables')
        values = {var: value for var, value in zip(vars, values)}
        res = Poly()
        for monomial, coeff in self.monomials.items():
            res += eval_monomial(monomial, coeff, values)
        return res


def eval_monomial(monomial, coeff, values):
    res = defaultdict(Perm)
    for var, exp in monomial.items():
        if var in values:
            coeff *= values[var]**exp
        else:
            res[var] = exp
    return Poly(defaultdict(Perm, {frozendict(res): coeff}))


def Var(name):
    return Poly(defaultdict(Perm, {frozendict({name: 1}): 1}))


X = Var('X')
Y = Var('Y')
Z = Var('Z')


from itertools import count
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint


def minimal_poly(elem):
    if elem == 0:
        return Poly(1)
    if isinstance(elem, int):
        elem = elem * Perm(1)
    for deg in count(0):
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
        c = np.zeros(deg + 1)
        constraints = LinearConstraint(A, L, U)
        integrality = np.ones(deg + 1)
        res = milp(c=c, bounds=Bounds(),
                   constraints=constraints,
                   integrality=integrality)
        if res.x is not None:
            return sum(int(n) * X**i for n, i in zip(res.x, count(0)))


from sympy import Matrix

def minimal_poly_nullspace(elem):
    if isinstance(elem, int):
        elem = elem * Perm(1)
    for deg in count(0):
        lengths_values = (elem**deg).lengths()
        # print(lengths_values)
        rows = len(lengths_values)
        # print(lengths_values)
        A = np.zeros((rows, deg+1),dtype=int)
        for j in range(0,deg+1):
            potenza_j = elem**j
            for i in range(rows):
                # print(elem**j,(elem**j).coeff(lengths_values[i]))
                A[i][j] = (potenza_j).coeff(lengths_values[i])

        B = Matrix(A)
        #print(B)
        
        solspace = B.nullspace()

        # print(solspace)
        pols = [Poly(*[int(n) for n in reversed(sol)]) for sol in solspace if sol[-1]!=0 ]
        if pols != []:
            return pols


def linear_combination(elem, deg):
    rows = (elem**deg).lengths()[-1] +1 
    A = np.zeros((rows, deg))
    for i in range(rows):
        for j in range(deg):
            A[i][j] = (elem**j).coeff(i)
    
    L = np.zeros(rows)
    U = np.zeros(rows)
    for i in range(rows):
       L[i] = (elem**deg).coeff(i)
       U[i] = (elem**deg).coeff(i)
    
    constraints = LinearConstraint(A, L, U)
    
    # print(A,L,U)

    c = np.zeros(deg)
    integrality = np.ones(deg)
    
    res = milp(c=c, bounds=Bounds(lb=0),
               constraints=constraints,
               integrality=integrality)
    # print(res)
    if res.x is None:
        # print(res)
        return None
    return Poly(*([-1]+[int(n) for n in reversed(res.x)]))
