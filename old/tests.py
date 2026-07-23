from permutations import *

# Minimal poly of C(2) + C(3) + C(5)

P_2_3_5 = (np.array([1, -35, 487, -3425, 12736, -23540, 16800, 0]) @
           [X**i for i in reversed(range(8))])

# Minimal poly of C(2) + C(3) + C(7)

P_2_3_7 = (np.array([1, -48, 946, -9864, 58345, -194136, 333540, -226800, 0]) @
           [X**i for i in reversed(range(9))])

# Minimal poly of 1 + C(3) (?)

P_1_3 = X**2 - 5*X + 4

# Tests

x0 = variable('x0')
x1 = variable('x1')
x2 = variable('x2')
x3 = variable('x3')
x4 = variable('x4')
A = C(2)+C(3)
PA = x4*A**4 + x3*A**3 + x2*A**2 + x1*A + x0


# Testing if all minimal polynomials of permutations split over N

def is_int(A):
    return isinstance(A, int) or A.cycles() == {C(1)}


def test_splitting(A):
    E = A.minimal_equation()
    sols = E.solutions_univariate()
    int_sols = [sol for sol in sols if is_int(sol)]
    return len(int_sols) == E.degree()


def test_all_splittings(size=None):
    if not size:
        sizes = count(1)
    else:
        sizes = [size]
    for size in sizes:
        print(size)
        for A in Permutation.generate(size):
            if not test_splitting(A):
                print(A.minimal_equation())
            assert test_splitting(A), f'{A} failed'



from math import ceil


# k-permutation = C(n[1]) + ... + C(n[k]) with n[i] != n[j] for i != j

def k_permutations(size, k):
    for A in Permutation.generate(size):
        if (len(A.cycles()) == k and
            all(A.multiplicity(length) == 1
                for length in A.lengths())):
            yield A


# k-permutations without 1
            
def good_k_permutations(size, k):
    for A in k_permutations(size, k):
        if 1 not in A.lengths():
            yield A


def pairwise_coprime(nums):
    return all(gcd(a, b) == 1
               for (a, b) in combinations(nums, 2))


def expected_degree(A):
    return 2**len([c for c in A.cycles() if c != 1])


def test_k_permutations(k):
    for size in count(1):
        for A in good_k_permutations(size, k):
            if pairwise_coprime(A.lengths()):
                degree = A.minimal_equation().degree()
                if expected_degree(A) != degree:
                    print(A, '->', degree)

# Apparently the degree of C(2) + C(3) + ... + C(p[n]) is https://oeis.org/A082548
# Indeed the roots seem to be {sum(X) for X in powerset(cycles of the equation)}


def sums_of_subsets(S):
    return {sum(T) for T in powerset(S)}


def equation_with_roots(R):
    P = prod((X-n) for n in R)
    return Equation(P, 0)


# Check if the minimal polynomial of C(2) + C(3) + C(5) + ... + C(p_n)
# has integer roots exactly sums_of_subsets({2, 3, 5, ..., p_n}): true
# until C(2) + ... + C(23)


from math import isqrt


def is_prime(n):
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def primes():
    for n in count(2):
        if is_prime(n):
            yield n


def test_min_poly_roots_primes():
    A = 0
    N = set()
    for p in primes():
        A += C(p)
        N.add(p)
        E = A.minimal_equation()
        S = sums_of_subsets(N)
        for s in S:
            assert E.P(s) == E.Q(s)
        assert len(S) == E.degree()
        print(A)


# Check if the minimal polynomial of C(n_1) + ... + C(n_m) with n_i
# pairwise coprime has integer roots exactly sums_of_subsets({n_1, ..., n_m})
# true for n_1 + ... + n_m <= 67


def coprime_partitions(n):
    for p in partitions(n):
        if 1 not in p and pairwise_coprime(p):
            yield p


def test_min_poly_roots_coprimes():
    for n in count(1):
        for part in coprime_partitions(n):
            N = set(part)
            A = Permutation.of(sum(C(k) for k in part))
            E = A.minimal_equation()
            S = sums_of_subsets(N)
            var = E.vars()[0]
            for s in S:
                assert E.P(**{var: s}) == E.Q(**{var: s})
            assert len(S) == E.degree()
            print(n, A)


# Check if the minimal polynomial of C(n_1) + ... + C(n_m)
# without hypotheses on n_1, ..., n_m
# has integer roots exactly sums_of_lcm_dict({n_1, ..., n_m})
# true for n_1 + ... + n_m <= 29


def lcm_dict(S):
    d = {}
    for T in powerset(S):
        l = lcm(*T)
        # d[l] = d.get(l, []) + [set(T)]
        # the subsets of S are generated in order of <=
        # the union of sets with the same lcm has the same lcs
        # so we only take the largest one
        d[l] = T
    return d


def sums_of_lcm_dict(S):
    return {sum(T) for T in lcm_dict(S).values()}


def test_min_poly_roots():
    for n in count(1):
        for part in partitions(n):
            A = Permutation.of(sum(C(k) for k in part))
            N = part
            E = A.minimal_equation()
            S = sums_of_lcm_dict(N)
            var = E.vars()[0]
            for s in S:
                assert E.P(**{var: s}) == E.Q(**{var: s})
            assert len(S) == E.degree()
            print(f'size {n}\tdeg = {E.degree()}\tncycles = {len(A.cycles())}\t\t{A}')


def test_min_poly_degrees():
    degs_max = []
    for n in count(1):
        deg_max = 0
        for part in partitions(n):
            A = Permutation.of(sum(C(k) for k in part))
            E = A.minimal_equation()
            deg = E.degree()
            # S = sums_of_lcm_dict(part)
            # deg = len(S)
            deg_max = max(deg_max, deg)
        if deg_max not in degs_max:
            degs_max.append(deg_max)
        print(n, degs_max)


def test_min_poly_degrees_by_ncycles():
    d = {}
    for n in count(1):
        for part in partitions(n):
            A = Permutation.of(sum(C(k) for k in part))
            # N = part
            E = A.minimal_equation()
            # S = sums_of_lcm_dict(N)
            # var = E.vars()[0]
            # for s in S:
            #     assert E.P(**{var: s}) == E.Q(**{var: s})
            # assert len(S) == E.degree()
            ncycles = len(A.cycles())
            deg = E.degree()
            if d.get(ncycles, np.inf) > deg:
                print(A, '->', ncycles)
            d[ncycles] = min(d.get(ncycles, np.inf), deg)
        print(f'{n}: {d}')


# def equation_from_polynomial(R):
#     P = Polynomial.of(0)
#     Q = Polynomial.of(0)
#     for mono, coeff in R.terms().items():
#         coeff = Permutation.of(coeff)
#         mult = list(coeff._cycles.values())[0]
#         if mult > 0:
#             P += coeff * Polynomial.of(mono)
#         else:
#             Q -= coeff * Polynomial.of(mono)
#     return Equation(P, Q)


def equation_from_polynomial(R):
    P = Polynomial.of(0)
    Q = Polynomial.of(0)
    for mono, coeff in R.terms().items():
        coeff = Permutation.of(coeff)
        for length, mult in coeff._cycles.items():
            term = mult * C(length) * Polynomial.of(mono)
            if mult > 0:
                P += term
            else:
                Q -= term
    return Equation(P, Q)


# Tests

A = C(2) + C(3) + C(5) + C(7)


# Counterexample to glueing by multiplying by a cycle

E1 = Equation(X + Y, C(3))
E2 = Equation(C(2) + Z, C(2) * X)

wrongE = Equation(C(5) * E1.P + C(7) * E2.P,
                  C(5) * E1.Q + C(7) * E2.Q)

E = equation_from_polynomial((C(5) - 5) * (E1.P - E1.Q) +
                             (C(7) - 7) * (E2.P - E2.Q))


# Test that P(X) + Q(Y) == 0 if and only if P(X) == Q(Y) == 0
# with P, Q ∈ P[-1][X] but X, Y ∈ P

def test_poly_sum_zero(P, Q, bound):
    E = equation_from_polynomial(P + Q)
    print(E)
    for sizeA in range(bound + 1):
        print(sizeA)
        for A in Permutation.generate(sizeA):
            for sizeB in range(bound + 1):
                for B in Permutation.generate(sizeB):
                    if E.eval(A, B) and P(A) != 0:
                        print('FALSE!', A, B)


# Test that P(X) + Q(Y) == 0 if and only if P(X) == Q(Y) == 0
# with P, Q ∈ P[-1][X] AND X, Y ∈ P[-1]

def test_poly_sum_zero_extended(P, Q, bound):
    E = equation_from_polynomial(P + Q)
    print(E)
    for n1, n2, n3, n4 in product(range(bound + 1), repeat=4):
        print(n1, n2, n3, n4)
        for A1, A2, A3, A4 in product(Permutation.generate(n1),
                                      Permutation.generate(n2),
                                      Permutation.generate(n3),
                                      Permutation.generate(n4)):
            if E.eval(A1 - A2, A3 - A4) and P(A1 - A2) != 0:
                print(f'FALSE! {A1}, {A2}, {A3}, {A4} -> {A1-A2}, {A3-A4}')


# Boolean formula (x2 v x3 v x5) & (~x3 v ~x5 v ~x7)

A = C(2) + C(3) + C(5) + C(7)
# X <= A
# C(2*3*5) <= C(2*3*5) * X
# C(3*5*7)*(X + 1) <= C(3*5*7)*A

# E = equation_from_polynomial((C(11) - 11) * (X - A) +
#                              (C(13) - 13) * (C(2*3*5) - C(2*3*5) * X) +
#                              (C(17) - 17) * (C(3*5*7)*(X + 1) - C(3*5*7)*A) +
#                              Y)


# Natural system:
# Y + 1   == X
# Y + 2*X == 5
# has solution (2, 1)

E = equation_from_polynomial((C(2) - 2) * (Y + 1 - X) +
                             (C(3) - 3) * (Y + 2*X - 5))
# E = equation_from_polynomial((C(3) - 2) * (Y + 1 - X) +
#                              (3 + (-1)*C(2)) * (Y + 2*X - 5))


# Natural system:
# Y + 4*X == X**2 + 6
# Y = X + 2
# has solutions (1, 3) and (4, 6)

E1 = equation_from_polynomial((C(5) - 5) * (Y + 4*X - X**2 - 6) +
                              (C(11) - 11) * (Y - X - 2))


def square(side):
    return product(range(side + 1), repeat=2)


def test_equation(E, domain):
    for a, b in domain:
        print(f'{a=}, {b=}')
        for A, B in product(Permutation.generate(a),
                            Permutation.generate(b)):
            if E.eval(**{'X': A, 'Y': B}):
                print(f'SOLUTION! X={A}, Y={B}')


# Natural system:
# Y + 4*X == X**2 + 6
# Y = X + 1
# has no solution

E2 = equation_from_polynomial((C(2) - 2) * (Y + 4*X - X**2 - 6) +
                              (C(3) - 3) * (Y - X - 1))


# Natural system:
# Y + X**2 + 4*X = 16
# 2*Y = 3*X + 2
# has one solution (2, 4)

E3 = equation_from_polynomial((C(2) - 2) * (Y + X**2 + 4*X - 16) +
                              (C(7) - 7) * (2*Y - 3*X - 2))

# Natural system:
# Y + X**2 + 4*X = 16
# 2*Y = 3*X + 3
# has no solution

E4 = equation_from_polynomial((C(2) - 2) * (Y + X**2 + 4*X - 16) +
                              (C(7) - 7) * (2*Y - 3*X - 3))

E5 = equation_from_polynomial((C(13) - 13) * (Y + X**2 + 4*X - 16) +
                              (C(17) - 17) * (2*Y - 3*X - 3))


# Boolean formula (x2 v x3 v x5) & (~x3 v ~x5 v ~x7)

A = C(2) + C(3) + C(5) + C(7)
# X <= A
# C(2*3*5) <= C(2*3*5)*X
# C(3*5*7)*(X + 1) <= C(3*5*7)*A

E6 = equation_from_polynomial((X - A)**2 +
                              (C(2*3*5) - C(2*3*5)*X)**2 +
                              (C(3*5*7)*(X + 1) - C(3*5*7)*A)**2 +
                              Y)


# Equations with glueing with (C(3)-2) and (3-C(2))

G1 = C(3) - 2
G2 = 3 + (-1)*C(2)
# G1 = C(7) - 3
# G2 = 7 + (-1)*C(3)

E7 = equation_from_polynomial(G1 * (Y + 1 - X) +
                              G2 * (Y + 2*X - 6))
E8 = equation_from_polynomial(G1 * (Y + 1 - X) +
                              G2 * (Y + 2*X - 5))
E9 = equation_from_polynomial(G1 * (Y + 1 - X) +
                              G2 * (Y + 2*X + Z - 5))
E10 = equation_from_polynomial(G1 * (Y + 4*X - X**2 - 6) +
                               G2 * (Y - X - 1))
E11 = equation_from_polynomial(G1 * (Y + X**2 + 4*X - 16) +
                               G2 * (2*Y - 3*X - 2))
E12 = equation_from_polynomial(G1 * (X - Y) +
                               G2 * 2)


# Simulating matrices

a = variable('a')
b = variable('b')
c = variable('c')
d = variable('d')
M = a*C(2) + b*C(4) + c*C(6) + d*C(8)


# Irreducibility test

def is_irreducible(A):
    A = Permutation.of(A)
    if A == 0 or A == 1:
        return False
    n = len(A)
    found = False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            m = n // d
            for B1 in Permutation.generate(d):
                for B2 in Permutation.generate(m):
                    if B1 * B2 == A:
                        print(f'{B1}, {B2}')
                        found = True
    return not found
