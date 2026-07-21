class Permutations(CombinatorialFreeModule):

    def __init__(self):
        CombinatorialFreeModule.__init__(
            self, ZZ, NN, prefix='C',
            category=Category.join(
                (AlgebrasWithBasis(ZZ),
                 CommutativeRings())))

    def product_on_basis(self, left, right):
        return gcd(left, right) * C[lcm(left, right)]

    def one_basis(self):
        return 1

    def _repr_(self):
        return '(Semi)ring of permutations'


PP = Permutations()
C = PP.basis()
R.<X> = PP[]
