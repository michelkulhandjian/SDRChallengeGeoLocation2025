import numpy as np
from numpy.fft import fft, ifft, fftshift
c = 299792458.0

def xcorr_fft(a,b, up=8):
    N = len(a)+len(b)-1
    M = 1<<int(np.ceil(np.log2(N*up)))
    A = fft(a, M); B = fft(b, M)
    R = A*np.conj(B)
    r = ifft(R)
    r = np.concatenate([r[-(M//2):], r[:M//2]])
    return np.abs(r), np.linspace(-(M//2), M//2-1, M)/up

def toa_est(iq, preamble):
    r, lags = xcorr_fft(iq, preamble, up=8)
    k = np.argmax(r)
    # quadratic interp
    if 1 <= k < len(r)-1:
        y0, y1, y2 = r[k-1], r[k], r[k+1]
        denom = (y0 - 2*y1 + y2)
        delta = 0.5*(y0 - y2)/denom if denom != 0 else 0.0
    else:
        delta = 0.0
    return (k + delta)  # in upsampled samples

def tdoa_solve_2d(anchors, tdoas, ref, fs_up):
    # Nonlinear LS (Gauss-Newton) for 2D, anchors Nx2, tdoas relative to ref
    x = np.mean(anchors, axis=0)  # init at centroid
    for _ in range(15):
        H = []; r = []
        for i,(xi,yi) in enumerate(anchors):
            if i==ref: continue
            di = np.linalg.norm(x - anchors[i])
            dr = np.linalg.norm(x - anchors[ref])
            g = di - dr
            tau = tdoas[i]  # seconds relative to ref (0 for ref)
            r.append(g - c*tau)
            # Jacobian
            gi = (x - anchors[i])/(di+1e-9) - (x - anchors[ref])/(dr+1e-9)
            H.append(gi)
        H = np.vstack(H); r = np.array(r)
        dx, *_ = np.linalg.lstsq(H, -r, rcond=None)
        x = x + dx
        if np.linalg.norm(dx) < 1e-3:
            break
    return x

# Example usage:
FS = 25e6
UPS = 8
bias = np.array([0.0, 12.5e-9, -7e-9, 5e-9])  # per-channel bias (seconds) from calibration
anchors = np.array([[0,0],[15,0],[15,12],[0,12]])  # meters
# Load IQ per channel, estimate TOA in samples_upsampled:
# toa_samp[i] = toa_est(iq_i, preamble)
# Convert to seconds: (toa_samp / (FS*UPS)) - bias[i]
# Make tdoas rel. to ref channel 0:
# tdoa[i] = toa_sec[i] - toa_sec[0]
# Then solve:
# pos = tdoa_solve_2d(anchors, tdoa, ref=0, fs_up=FS*UPS)
