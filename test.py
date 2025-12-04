from NTRUSolve import *
from fft import *

p = gen_poly(8)

print(f"Inicio: poli = {p}, len = 8\n\n")

fft_p = fft(p)

print(f"\n\nfft: poli = {fft_p}, len = {len(fft_p)}\n\n")

pp = ifft(fft_p)

print(f"\n\nifft: poli = {pp}, len = {len(pp)}\n\n")

print(adj(p))
print(adj_fft(fft_p))