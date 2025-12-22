from common import split, merge         
from fft_constants import roots_dict 

# separa um polinomio f na fft em 2
def split_fft(f_fft):
    n = len(f_fft)
    w = roots_dict[n]
    f0_fft = [0] * (n // 2)
    f1_fft = [0] * (n // 2)
    for i in range(n // 2):
        f0_fft[i] = 0.5 * (f_fft[2 * i] + f_fft[2 * i + 1])
        f1_fft[i] = 0.5 * (f_fft[2 * i] - f_fft[2 * i + 1]) * w[2 * i].conjugate()
    return [f0_fft, f1_fft]

# agrupa dois polinômios na fft
def merge_fft(f_list_fft):
    f0_fft, f1_fft = f_list_fft
    n = 2 * len(f0_fft)
    w = roots_dict[n]
    f_fft = [0] * n

    for i in range(n // 2):
        f_fft[2 * i + 0] = f0_fft[i] + w[2 * i] * f1_fft[i]
        f_fft[2 * i + 1] = f0_fft[i] - w[2 * i] * f1_fft[i]

    return f_fft

# calcula a fft de um polinomio mod (x^n + 1)
def fft(f):
    #print(f"poli = {f}, len = {len(f)}")
    n = len(f)

    if n > 2:
        f0, f1 = split(f)
        f0_fft = fft(f0)
        f1_fft = fft(f1)
        f_fft = merge_fft([f0_fft, f1_fft])

    elif n == 2:
        f_fft = [0] * n
        f_fft[0] = f[0] + 1j * f[1]
        f_fft[1] = f[0] - 1j * f[1]

    return f_fft

# calcula a fft inversa de um polinomio mod (x^n + 1)
def ifft(f_fft):
    #print(f"i poli = {f_fft}, len = {len(f_fft)}")
    n = len(f_fft)

    if (n > 2):
        f0_fft, f1_fft = split_fft(f_fft)
        f0 = ifft(f0_fft)
        f1 = ifft(f1_fft)
        f = merge([f0, f1])

    elif (n == 2):
        f = [0] * n
        f[0] = f_fft[0].real
        f[1] = f_fft[0].imag

    return f

# multiplca dois polinomios na representação fft
def mul_fft(f_fft, g_fft):
    deg = len(f_fft)
    return [f_fft[i] * g_fft[i] for i in range(deg)]

# divite dois polinomios na representação fft
def div_fft(f_fft, g_fft):
    assert len(f_fft) == len(g_fft)
    deg = len(f_fft)
    return [f_fft[i] / g_fft[i] for i in range(deg)]

# multiplica dois polinomios
def mul(f, g):
    return ifft(mul_fft(fft(f), fft(g)))

# divide dois polinomios
def div(f, g):
    return ifft(div_fft(fft(f), fft(g)))

# soma dois polinômios
def add(f, g):
    assert len(f) == len(g)
    n = len(f)
    return [f[i] + g[i] for i in range(n)]

# subtrai dois polinômios
def sub(f, g):
    assert len(f) == len(g)
    n = len(f)
    return [f[i] - g[i] for i in range(n)]

# 
def adj_fft(f_fft):
    """Ajoint of a polynomial (FFT representation)."""
    deg = len(f_fft)
    return [f_fft[i].conjugate() for i in range(deg)]

# soma dois polinomios
def adj(f):
    return ifft(adj_fft(fft(f)))

# soma dois polinomios na fft
def add_fft(f_fft, g_fft):
    return add(f_fft, g_fft)

# subtrai dois polinomios na fft
def sub_fft(f_fft, g_fft):
    return sub(f_fft, g_fft)

# multiplica um polinomio por um escalar k modulo m
def scalar_mul(f, k, m):
    return [(k*i) % m for i in f]

def poly_mod(a, m):
    return [x % m for x in a]

# razão entre:
#   * O grau de n
#   * o número de coeficientes complexos da ntt
fft_ratio = 1