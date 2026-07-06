from itertools import product, chain

# E = n | X(n) | C(n) | E + E | E * E | E**n


class Expr:
    pass


class C(Expr):

    def __init__(self, n):
        self.n = n

    def __repr__(self):
        return f'C({self.n})'

    def eval(self, ass={}):
        return Plus(self)

    def __mul__(self, other):
        return Times(self, other)

    def __rmul__(self, other):
        return Times(other, self)

    def __add__(self, other):
        return Plus(self, other)

    def __radd__(self, other):
        return Plus(other, self)


class X(Expr):

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f'X({self.name})'

    def eval(self, ass={}):
        if self.name in ass:
            return Plus(ass[self.name])
        else:
            return Plus(self)


class Plus(Expr):

    def __init__(self, *terms):
        self.terms = terms

    def __repr__(self):
        return f'Plus{self.terms}'

    def eval(self, ass={}):
        terms = (term.eval(ass) for term in self.terms)
        return Plus(*chain(*(term.terms for term in terms)))


class Times(Expr):

    def __init__(self, *terms):
        self.terms = terms

    def __repr__(self):
        return f'Times{self.terms}'

    def eval(self, ass={}):
        terms = (term.eval(ass) for term in self.terms)
        prod = product(*(term.terms for term in terms))
        return Plus(*(Times(*p) for p in prod))


E = Times(Plus(X(2), C(2)),
          Plus(C(3), X(1)))
