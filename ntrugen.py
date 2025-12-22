from common import sqnorm
from samplerz import samplerz
from fft import *
from ntt import *

q = 12 * 1024 + 1

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

def field_norm(f):
    n2 = len(f) // 2
    fe = [f[2 * i] for i in range(n2)]
    fo = [f[2 * i + 1] for i in range(n2)]
    fe2 = mult_mod(fe, fe)
    fo2 = mult_mod(fo, fo)
    res = fe2[:]
    
    for i in range(n2 - 1):
        res[i + 1] += fo2[i]
    res[0] += fo2[n2 - 1]
    return res

# MDC extendido de f0 e g0
def mdc_x(b, n):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return b, x0, y0

def bitsize(a):
    val = abs(a)
    res = 0
    while val:
        res += 8
        val >>= 8
    return res

# Reduz (F, G) em relação a (f, g)
# Feito usando redução de Babai
def reduce(f, g, F, G):
    n = len(f)
    size = max(53, bitsize(min(f)), bitsize(max(f)), bitsize(min(g)), bitsize(max(g)))

    f_adjust = [elt >> (size - 53) for elt in f]
    g_adjust = [elt >> (size - 53) for elt in g]
    fa_fft = fft(f_adjust)
    ga_fft = fft(g_adjust)

    while(1):
        Size = max(53, bitsize(min(F)), bitsize(max(F)), bitsize(min(G)), bitsize(max(G)))
        if Size < size:
            break

        F_adjust = [elt >> (Size - 53) for elt in F]
        G_adjust = [elt >> (Size - 53) for elt in G]
        Fa_fft = fft(F_adjust)
        Ga_fft = fft(G_adjust)

        den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)), mul_fft(ga_fft, adj_fft(ga_fft)))
        num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)), mul_fft(Ga_fft, adj_fft(ga_fft)))
        k_fft = div_fft(num_fft, den_fft)
        k = ifft(k_fft)
        k = [int(round(elt)) for elt in k]
        if all(elt == 0 for elt in k):
            break
        fk = mult_mod(f, k)
        gk = mult_mod(g, k)
        for i in range(n):
            F[i] -= fk[i] << (Size - size)
            G[i] -= gk[i] << (Size - size)
    return F, G

def gen_poly(N):
    #1.17 * sqrt(12289 / 8192)
    sigma = 1.43300980527773

    assert(N < 4096)
    f0 = [samplerz(0, sigma, sigma - 0.001) for _ in range(4096)]
    f = [0] * N
    k = 4096 // N
    for i in range(N):
        # We use the fact that adding k Gaussian samples of std. dev. sigma
        # gives a Gaussian sample of std. dev. sqrt(k) * sigma.
        f[i] = sum(f0[i * k + j] for j in range(k))
    return f

def ntru_solve(f, g):
    print(f"F = {f}, G = {g}")
    n = len(f)
    if n == 1:
        f0 = f[0]
        g0 = g[0]
        d, u, v = mdc_x(f0, g0)
        if d != 1:
            raise ValueError
        else:
            return [-q * v],[q * u]
    else:
        fp = field_norm(f)
        gp = field_norm(g)
        Fp, Gp = ntru_solve(fp, gp)
        F = mult_mod(lift(Fp), galois_conjugate(g))
        G = mult_mod(lift(Gp), galois_conjugate(f))
        F, G = reduce(f, g, F, G)
        return F, G

# Calcula a norma quadrada de Gram-Schmidt da matriz do NTRU gerada por f e g
def gs_norm(f, g, q):
    sqnorm_fg = sqnorm([f, g])
    ffgg = add(mul(f, adj(f)), mul(g, adj(g)))
    Ft = div(adj(g), ffgg)
    Gt = div(adj(f), ffgg)
    sqnorm_FG = (q ** 2) * sqnorm([Ft, Gt])
    return max(sqnorm_fg, sqnorm_FG)

# Gera os polinomios f, g, F e G satisfazendo a equação do NTRU
def ntru_gen(n):
    while True:
        f = gen_poly(n)
        g = gen_poly(n)
        if gs_norm(f, g, q) > (1.17 ** 2) * q:
            continue
        f_ntt = ntt(f)
        if any((elem == 0) for elem in f_ntt):
            continue
        try:
            F, G = ntru_solve(f, g)
            F = [int(coef) for coef in F]
            G = [int(coef) for coef in G]
            return f, g, F, G

        except ValueError:
            continue
    
def main():
    N = 4

    while True:
        f = gen_poly(N)
        g = gen_poly(N)

        print(f"\n\nInicio: f = {f}, g = {g}\n\n")

        try:
            F, G = ntru_solve(f, g)
            print(f"\n\nFinal: F = {F}, G = {G}")
            return
        
        except ValueError:
            continue

#main()