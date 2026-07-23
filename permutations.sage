from sage.functions.log import logb


class Permutation(CombinatorialFreeModule.Element):

    def cycles(self):
        return sorted(term.leading_term().leading_monomial()
                      for term in self.terms())

    def is_cycle(self):
        return len(self) == 1 and self.leading_coefficient() == 1

    def size(self):
        return int(sum(multiplicity * length
                       for length, multiplicity in self.items()))

    def is_irreducible(self):
        if self == 0 or self == 1:
            return False
        return all(A * B != self
                   for m, n in _proper_divisor_pairs(self.size())
                   for A in PP.of_size(m)
                   for B in PP.of_size(n))

    def sqrt(self, n=2):
        root = PP(0)
        power = PP(0)
        while power.size() < self.size():
            root += min((self - power).cycles())
            power = root^n
        if power != self:
            return None
        return root

    def factor(self):
        if self == 0:
            raise ArithmeticError('factorization of 0 is not defined')
        if self == 1:
            return Factorization([])
        if self.is_cycle():
            # Optimisation (notably for use with divisors())
            F = factor(self.size())
            return Factorization([(C[b^e], 1) for b, e in F])
        for m, n in _proper_divisor_pairs(self.size()):
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
    def up_closure(generators):
        return sorted(set(PP.prod(S).leading_monomial()
                          for S in powerset(generators)))

    @staticmethod
    def down_closure(generators):
        return sorted(set(div for gen in generators
                          for term in gen.terms()
                          for cycle in term.cycles()
                          for div in divisors(cycle)))

    # @staticmethod
    # def solve_linear(P):
    #     if P.degree() != 1:
    #         raise ValueError(f'{P} is not linear')
    #     D = Permutations.down_closure(P.coefficients())
    #     B = Permutations.up_closure(D)
    #     return B
    #     # TODO here

    @staticmethod
    def solve_univariate(P, all=False):
        if len(P.variables()) > 1:
            raise ValueError(f'{P} is not univariate')
        if (len(P.terms()) == 2 and min(P.exponents()) == 0
              or len(P.terms()) == 1) and P.leading_coefficient() == 1:
            # P == X^n - A
            A = -P.constant_coefficient()
            return A.sqrt(P.degree())
        identity = lambda i: i
        cardinality = PP.module_morphism(identity, codomain=ZZ)
        q = P.map_coefficients(cardinality)
        roots = q.roots(multiplicities=False)
        solutions = (A for size in roots
                     for A in PP.of_size(size)
                     if P(A) == 0)
        if all:
            return list(solutions)
        solution = next(solutions, None)
        if solution:
            return [solution]
        return []

    @staticmethod
    def solve_pseudo_injective(P, all=False):
        if all:
            raise NotImplementedError(
                'enumeration of all solutions not implemented yet')
        B = -P.constant_coefficient()
        P += B
        if not _is_pseudo_injective(P):
            raise ValueError(f'{P} is not pseudo-injective')
        X = PP(0)
        s = _seed(P)
        while P(X) != B:
            X += C[_alcm(s, B - P(X))]
            if P(X).size() > B.size() or not P(X) <= B:
                return []
        return [X]


PP = Permutations()

C = PP.basis()

_R.<X> = PP[]


def _proper_divisor_pairs(n):
    return ((d, n // d)
             for d in range(2, isqrt(n) + 1)
             if n % d == 0)


def _alcm(a, b):
    if a not in NN:
        a = min(a.cycles()).size()
    if b not in NN:
        b = min(b.cycles()).size()
    if b % a != 0:
        raise ValueError(f'{a} does not divide {b}')
    k = ceil(logb(9, 2))
    return gcd(NN(b/a)^k, b)


def _cycles(P):
    return sorted(x for c in P.coefficients()
                  for x in c.cycles())


def _is_pseudo_injective(P):
    cycles = _cycles(P)
    length = cycles[0].size()
    return all(cycle.size() % length == 0
               for cycle in cycles)


def _seed(P):
    cycles = _cycles(P)
    return cycles[0].size()


# Tests

_R.<Y, Z> = PP[]

P = C[3]*Y + (C[2] + C[7])*Z - C[6]
P = C[2]*X^2 + (C[4] + C[6])*X - 16*C[2] - 4*C[4] - 18*C[6] - C[12]
n = prod(primes_first_n(8))
P = (X^5 - C[n]^5) + (X - C[n])
