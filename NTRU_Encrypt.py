from ntrugen import *
from common import *
from math_utils import *
import math

class NTRUEncrypt:
    N = None
    p = None
    q = None

    f_poly = None
    g_poly = None
    h_poly = None

    f_p_poly = None
    f_q_poly = None

    R_poly = None

    def __init__(self, N, p, q):
        self.N = N
        self.p = p
        self.q = q

    def gen_poly(self, d, neg = 0):
        zeros = [0] * (self.N - 2 * d - neg)
        ones = [1] * d
        minus_ones = [-1] * (d + neg)

        poly = zeros + ones + minus_ones
        return csprng_shuffle(poly)
    
    def invert_poly(self, f, x):
        if is_2_power(x):
            return invert_poly_2power(f, g)
        else:
            raise NotImplementedError("Somente q = 2^k implementado")

    def generate_random_keys(self):
        g_poly = self.gen_poly(int(math.sqrt(self.q)))
        tries = 1000000000

        while tries > 0 and self.h_poly is None:
            f_poly = self.gen_poly(self.N // 3, neg=-1)

            try:
                self.generate_public_key(f_poly, g_poly)
            except Exception as e:
                tries -= 1

        if self.h_poly is None:
            raise Exception("Couldn't generate invertible f")

    def generate_public_key(self, f_poly, g_poly):
        self.f_poly = f_poly
        self.g_poly = g_poly

        self.f_p_poly = self.invert_poly(self.f_poly, self.p)    
        self.f_q_poly = self.invert_poly(self.f_poly, self.q)

        p_f_q_poly = scalar_mul(self.f_q_poly, self.p, self.q)

        self.h_poly = poly_mod(mult_mod(p_f_q_poly, g_poly), self.q)

x = NTRUEncrypt(8, 3, 128)
x.generate_random_keys()

print(x.f_poly, x.g_poly, x.h_poly)
