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
