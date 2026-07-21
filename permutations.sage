class Permutation(CombinatorialFreeModule.Element):

    def cycles(self):
        return sorted(term.leading_term().leading_monomial()
                      for term in self.terms())

    def is_cycle(self):
        return len(self) == 1 and self.leading_coefficient() == 1

    def size(self):
        return sum(multiplicity * length
                   for length, multiplicity in self.items())

    def is_irreducible(self):
        if self == 0 or self == 1:
            return False
        size = self.size()
        for m, n in proper_divisor_pairs(size):
            for A in PP.of_size(m):
                for B in PP.of_size(n):
                    if A * B == self:
                        return False
        return True

    def factor(self):
        if self == 0:
            raise ArithmeticError(
                "factorization of 0 is not defined")
        if self == 1:
            return Factorization([])
        size = self.size()
        for m, n in proper_divisor_pairs(size):
            for A in PP.of_size(m):
                for B in PP.of_size(n):
                    if A * B == self:
                        return Factorization([(A, 1)]) * B.factor()
        return Factorization([(self,1)])


class Permutations(CombinatorialFreeModule):

    Element = Permutation

    def __init__(self):
        CombinatorialFreeModule.__init__(
            self, ZZ, NN, prefix='C',
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
    def of_size(size):
        for partition in Partitions(size):
            yield PP.sum(C[n] for n in partition)

    @staticmethod
    def irreducibles_of_size(size):
        for A in PP.of_size(size):
            if A.is_irreducible():
                yield A

    def __iter__(self):
        for size in NN:
            yield from self.of_size(size)


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
