import numpy as np
from numpy.fft import fft, ifft

def xcorr_fft(a, b, up=1):
    # cross-correlation of a (rx) with b (known), upsampled by 'up' via zero-padding
    N = len(a) + len(b) - 1
    M = 1 << int(np.ceil(np.log2(N)))
    A = fft(a, M); B = fft(b, M)
    R = A * np.conj(B)
    if up > 1:
        # zero-pad spectrum symmetrically
        pad = (up-1) * M
        Rzp = np.concatenate([R[:M//2], np.zeros(pad, dtype=R.dtype), R[M//2:]])
        r = ifft(Rzp)
        L = up * M
    else:
        r = ifft(R); L = M
    r = np.abs(r)
    # put peak near center (lag 0)
    r = np.concatenate([r[-(L//2):], r[:L//2]])
    return r

def quad_subsample(y, k):
    # 3-pt quadratic around index k
    if k <= 0 or k >= len(y)-1:
        return float(k), 0.0
    y1, y2, y3 = y[k-1], y[k], y[k+1]
    denom = (y1 - 2*y2 + y3)
    delta = 0.5*(y1 - y3)/denom if denom != 0 else 0.0
    return k + delta, delta

# usage:
# r = xcorr_fft(rx_block, known_preamble, up=16)
# k0 = np.argmax(r)
# k_hat, frac = quad_subsample(r, k0)
# tau_sec = (k_hat - len(r)//2) / (FS * up)   # convert lag index to seconds
