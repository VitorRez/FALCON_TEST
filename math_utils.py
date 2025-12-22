from common import *
from fft import *
import secrets

# pega um polinomio de tamanho n e transforma ele
# em um polinomio 2x maior com zeros entre cada coeficiente
def lift(f):
    n = len(f)
    res = [0] * (2 * n)
    for i in range(n):
        res[2 * i] = f[i]
    return res

# pega um polinomio e aplica uma transformação tal que
#   . termos de grau par ficam iguais
#   . termos de grau impar invertem o sinal
def galois_conjugate(f):
    n = len(f)
    return [((-1) ** i) * f[i] for i in range(n)]

# multiplicação entre dois polinômios sem fft
def mult(a, b):
    n = len(a)
    if n == 1:
        return [a[0] * b[0], 0]
    else:
        n2 = n // 2
        a0 = a[:n2]
        a1 = a[n2:]
        b0 = b[:n2]
        b1 = b[n2:]
        ax = [a0[i] + a1[i] for i in range(n2)]
        bx = [b0[i] + b1[i] for i in range(n2)]
        a0b0 = mult(a0, b0)
        a1b1 = mult(a1, b1)
        axbx = mult(ax, bx)
        for i in range(n):
            #print(f"i = {i} n = {n}")
            axbx[i] -= (a0b0[i] + a1b1[i])
        ab = [0] * (2 * n)
        for i in range(n):
            ab[i] += a0b0[i]
            ab[i + n] += a1b1[i]
            ab[i + n2] += axbx[i]
        return ab
    
# multiplicação entre polinomios seguido de uma redução mod(x^n + 1)
def mult_mod(a, b):
    n = len(a)
    ab = mult(a, b)
    abr = [ab[i] - ab[i + n] for i in range(n)]
    return abr

# MDC extendido de f0 e g0
def mdc_x(b, n):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return b, x0, y0

# 
def csprng_shuffle(p):
    p = p[:]
    for i in range(len(p) - 1, 0, -1):
        j = secrets.randbelow(i + 1)
        p[i], p[j] = p[j], p[i]
    return p

def is_prime(n):
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_2_power(n):
    return n != 0 and (n & (n - 1) == 0)

def invert_int(a, m):
    return pow(a, -1, m)

def poly_mod_int(f, m):
    return [x % m for x in f]

def poly_one(n):
    res = [0]*n
    res[0] = 1
    return res

def invert_poly_mod2(f):
    n = len(f)

    # caso base
    if n == 1:
        if f[0] % 2 == 0:
            raise ValueError("Polinômio não invertível em GF(2)")
        return [1]

    # split
    f0, f1 = split(f)

    # recursão
    g0 = invert_poly_mod2(f0)
    g1 = invert_poly_mod2(f1)

    # combina
    g0x = lift(g0)
    g1x = lift(g1)

    t0 = mult_mod(g0x, galois_conjugate(f))
    t1 = mult_mod(g1x, f)

    h = [(t0[i] + t1[i]) % 2 for i in range(n)]
    return h


def invert_poly_2power(f, q):
    n = len(f)

    # inversão inicial em GF(2)
    g = invert_poly_mod2([x % 2 for x in f])

    k = 1
    while (1 << k) < q:
        mod = 1 << (k + 1)

        fg = mult_mod(f, g)
        fg = poly_mod_int(fg, mod)

        two_minus_fg = sub(scalar_mul(poly_one(n), 2, mod), fg, mod)

        g = mult_mod(g, two_minus_fg)
        g = poly_mod_int(g, mod)

        k += 1

    return poly_mod_int(g, q)
