import numpy as np
from numpy.fft import fft, ifft

def xcorr_fft_upsampled(rx, ref, up=64):
    # FFT xcorr with frequency-domain zero-padding upsampling
    N = len(rx) + len(ref) - 1
    M = 1 << int(np.ceil(np.log2(N)))
    RX = fft(rx, M); REF = fft(ref, M)
    R = RX * np.conj(REF)
    if up > 1:
        pad = (up - 1) * M
        R = np.concatenate([R[:M//2], np.zeros(pad, dtype=R.dtype), R[M//2:]])
        L = up * M
    else:
        L = M
    r = ifft(R)
    r = np.abs(r)  # magnitude is stable for parabolic fit
    # center lags around zero
    r = np.concatenate([r[-(L//2):], r[:L//2]])
    return r, L

def quad_peak(y, k):
    # 3-point quadratic interpolation around index k
    if k <= 0 or k >= len(y)-1: return float(k)
    y1, y2, y3 = y[k-1], y[k], y[k+1]
    denom = (y1 - 2*y2 + y3)
    delta = 0.5 * (y1 - y3) / denom if denom != 0 else 0.0
    return k + delta

# usage:
# r, L = xcorr_fft_upsampled(rx_block, known_preamble, up=64)
# k0 = np.argmax(r)
# k_hat = quad_peak(r, k0)
# tau_sec = (k_hat - L//2) / (Fs * up)   # convert lag index to seconds
# -> subtract per-channel bias, form TDOAs, solve position.
