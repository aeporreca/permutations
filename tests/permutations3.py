from functools import total_ordering

# Frozendict should become part of Python 3.15
from frozendict import frozendict


@total_ordering
class Monomial:

    def __init__(self, vars={}):
        self.vars = frozendict(vars)

    @staticmethod
    def of(other):
        if isinstance(other, Monomial):
            return other
        if isinstance(other, str):
            return Monomial({other: 1})
        return Monomial({})

    def __hash__(self):
        return hash(self.vars)

    def __repr__(self):
        if not self.vars:
            return '1'
        vars_repr = []
        for var in sorted(self.vars.keys()):
            if self.vars[var] == 1:
                vars_repr.append(var)
            else:
                vars_repr.append(f'{var}**{self.vars[var]}')
        return ' * '.join(vars_repr)

    def __mul__(self, other):
        other = Monomial.of(other)
        vars = {}
        for var in self.vars.keys() | other.vars.keys():
            vars[var] = (self.vars.get(var, 0) +
                         other.vars.get(var, 0))
        return Monomial(vars)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        other = Monomial.of(other)
        return self.vars == other.vars

    def __le__(self, other):
        return repr(self) <= repr(other)


class Polynomial:

    def __init__(self, terms={}):
        self.terms = frozendict(terms)

    @staticmethod
    def of(coeff):
        if isinstance(coeff, Polynomial):
            return coeff
        if isinstance(coeff, Monomial):
            return Polynomial({coeff: 1})
        if coeff == 0:
            return Polynomial()
        return Polynomial({Monomial(): coeff})

    def __repr__(self):
        if not self.terms:
            return '0'
        terms_repr = []
        for monomial, coeff in self.terms.items():
            if monomial == 1:
                terms_repr.append(repr(coeff))
            elif coeff == 1:
                terms_repr.append(repr(monomial))
            else:
                terms_repr.append(f'{repr(coeff)} * {repr(monomial)}')
        return ' + '.join(terms_repr)

    def __add__(self, other):
        other = Polynomial.of(other)
        return Polynomial({monomial: self.terms.get(monomial, 0) +
                                     other.terms.get(monomial, 0)
                           for monomial in self.terms.keys() |
                                           other.terms.keys()})

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        other = Polynomial.of(other)
        return sum((self.terms.get(term1, 1) * self.terms.get(term2, 1) * term1 * term2
                    for term1 in self.terms
                    for term2 in other.terms),
                   start=Polynomial())

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exp):
        res = Polynomial.of(1)
        for i in range(exp):
            res *= self
        return res


def var(name):
    return Polynomial({Monomial({name: 1}): 1})


X = var('X')
Y = var('Y')
Z = var('Z')
