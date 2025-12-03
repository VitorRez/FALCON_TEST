q = 12 * 1024 + 1

# separa o polinomio f em 2 polinomios
def split(f):
    n = len(f)
    f0 = [f[2 * i] for i in range(n // 2)]
    f1 = [f[2 * i + 1] for i in range(n // 2)]
    return f0, f1

# agrupa dois polinômios em um só
def merge(f_list):
    f0, f1 = f_list
    n = 2 * len(f0)
    f = [0] * n
    for i in range(n // 2):
        f[2 * i] = f0[i]
        f[2 * i + 1] = f1[i]
    return f

# calcula a norma euclidiana quadrada de f
def sqnorm(f):
    norm = 0
    for i in f:
        for coef in i:
            norm += coef ** 2
    return norm