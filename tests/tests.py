from permutations import *

P = Poly(1, -2, 0)

# roots 0, 2, C₂, 3, C₃, 5, C₅, 2 + C₃, 3 + C₂, C₂ + C₃
Q = Poly(1, -10, 31, -30, 0)

# minimal_poly(C(4) + C(2), 3)
Poly(1, -8, 12, 0)

# minimal_poly(C(3) + C(9), 3)
Poly(1, -15, 36, 0)

# minimal_poly(C(2) + C(3), 4)
Poly(1, -10, 31, -30, 0)

def test(n):
    a = 2
    b = 3
    for i in range(n + 1):
        print(b**i - a**i)


def phi(a, b, n):
    return a**n - (b**n - a**n) * sum(a**i/(b**i - a**i)
               for i in range(2, n))


def test(r, n):
    return (r**n - 1) * (r**(n-1) - 1) / (r**(n-1) * (r - 1))


# def test(a, b, n):
#     return -b**(n-1) + a**(n-1) + (b**n - a**n) * sum(
#         ((b/a)**(k-1) - 1) / ((b/a)**(k-1) * b - a)
#         for k in range(2, n))

from math import floor

def test(a, b, n):
    return -b**n + a**n * (b - a + 1) + sum(
        (-a**k * (b - a + 1) + b**k) *
        floor((b**n - a**n) / (b**k - a**k))
        for k in range(2, n))


def P(X, Y, n):
    return ((1 + C(2))**n * X + (1 + C(3))**n * (X + 3 * Y)
            - (1 + C(2))**n * Y**2 - (1 + C(3))**n * 4)


P = Poly(1, 31, 0)
Q = Poly(0, 10 * (C(2) + C(3)), 30 * (C(2) + C(3)))
