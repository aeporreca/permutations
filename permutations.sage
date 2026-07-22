class Permutation(CombinatorialFreeModule.Element):

    def cycles(self):
        return sorted(term.leading_term().leading_monomial()
                      for term in self.terms())

    def is_cycle(self):
        return len(self) == 1 and self.leading_coefficient() == 1

    def size(self):
        return NN.sum(multiplicity * length
                      for length, multiplicity in self.items())

    def is_irreducible(self):
        if self == 0 or self == 1:
            return False
        return all(A * B != self
                   for m, n in proper_divisor_pairs(self.size())
                   for A in PP.of_size(m)
                   for B in PP.of_size(n))

    def factor(self):
        if self == 0:
            raise ArithmeticError(
                'factorization of 0 is not defined')
        if self == 1:
            return Factorization([])
        for m, n in proper_divisor_pairs(self.size()):
            for A in PP.irreducibles_of_size(m):
                for B in PP.of_size(n):
                    if A * B == self:
                        return Factorization([(A, 1)]) * B.factor()
        return Factorization([(self, 1)])


class Permutations(CombinatorialFreeModule):

    Element = Permutation

    def __init__(self):
        CombinatorialFreeModule.__init__(
            self, ZZ, PositiveIntegers(), prefix='C',
            category=Category.join(
                (AlgebrasWithBasis(ZZ),
                 CommutativeRings())))

    def product_on_basis(self, m, n):
        return gcd(m, n) * C[lcm(m, n)]

    def one_basis(self):
        return 1

    def _repr_(self):
        return '(Semi)ring of Permutations'

    @staticmethod
    def irreducibles():
        for size in NN:
            yield from PP.irreducibles_of_size(size)

    def __iter__(self):
        for size in NN:
            yield from PP.of_size(size)

    @staticmethod
    def of_size(size):
        for partition in Partitions(size):
            yield PP.sum(C[n] for n in partition)

    @staticmethod
    def irreducibles_of_size(size):
        for A in PP.of_size(size):
            if A.is_irreducible():
                yield A

    @staticmethod
    def solve_univariate(P):
        identity = lambda i: i
        f = PP.module_morphism(identity, codomain=ZZ)
        q = P.map_coefficients(f)
        roots = q.roots(multiplicities=False)
        return [A for size in roots
                for A in PP.of_size(size)
                if P(A) == 0]


PP = Permutations()

C = PP.basis()

R.<X> = PP[]


def proper_divisor_pairs(n):
    return ((d, n // d)
             for d in range(2, isqrt(n) + 1)
             if n % d == 0)


def up_closure(generators):
    return sorted(set(PP.prod(S).leading_monomial()
                      for S in powerset(generators)))
