from ntrugen import *

class NTRUEncrypt:

    def __init__(self, N):
        self.N = N

    def generate_random_keys(self):
        self.f, self.g, self.F, self.G = ntru_gen(self.N)
        self.h = div_zq(self.g, self.f)

    def encrypt()
    
x = NTRUEncrypt(1024)
x.generate_random_keys()
print(x.N, x.f, x.g, x.F, x.G, x.h)