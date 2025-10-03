#!/usr/bin/env python3
import time, json, argparse, math
import numpy as np
from numpy.fft import fft, ifft
from dataclasses import dataclass
from pathlib import Path

try:
    import uhd  # UHD Python bindings
except Exception as e:
    raise SystemExit("UHD Python bindings not found. Install UHD with Python or ensure 'uhd' is importable.") from e

# =======================
# ======= CONFIG =========
# =======================
Fs = 10e6                      # Sample rate (S/s)
FC = 2.45e9                    # Center frequency (Hz) - pick a quiet-ish channel
RX_GAIN = 30                   # dB (set to avoid clipping)
TX_GAIN = 0                    # dB baseband amplitude; RF gain controlled separately (set below)
TX_RF_GAIN = 10                # USRP TX gain in dB (adjust in your lab)
CHIP_RATE = 5e6                # 2 samples/chip at 10 MS/s
N_CHIPS = 10_000               # preamble length (chips)
UPS = 64                       # correlation upsample factor
CALIB_DURATION_S = 0.06        # seconds of capture (>= couple of preambles)
LOCATE_DURATION_S = 0.06
GUARD_SAMPS = int(0.5e-3 * Fs) # 0.5 ms guard between preambles
BIAS_FILE = "biases.json"
# Edit your device serials:
TX_SERIAL = None               # e.g., "30B0123" (None => pick first available)
RX_SERIALS = ["30A0001", "30A0002"]  # put your B210 serials here
# Receiver anchor positions (meters), same order as RX_SERIALS. Add as many as you have (>=3 for 2D fix)
ANCHORS = np.array([
    [0.0, 0.0],
    [15.0, 0.0],
    [15.0, 12.0],
    # add more rows if you have 4th receiver, etc.
], dtype=float)
# =======================

c = 299_792_458.0

@dataclass
class USRPHandle:
    u: any
    st: any

def generate_pn_bpsk(nchips:int, sps:int=2, seed:int=7):
    rng = np.random.default_rng(seed)
    chips = rng.integers(0, 2, size=nhips:=nchips).astype(np.int8)*2 - 1  # ±1
    # upsample (rectangular shaping)
    sig = np.repeat(chips.astype(np.float32), sps)
    return sig.astype(np.complex64)

def rrc_filter(x: np.ndarray, beta=0.25, span=8, sps=2):
    # Simple RRC shaping for demo (optional)
    # For brevity, rectangular is fine for indoor demo; skip if desired
    try:
        from scipy.signal import fir_filter_design as ffd  # noqa
    except Exception:
        return x  # fallback: rectangular if SciPy not present
    from scipy.signal import firwin, upfirdn
    N = span * sps * 2 + 1
    # crude RRC via root-raised-cosine freq shape is omitted for brevity; use rectangular or add proper rrc design
    return x  # keep rectangular here (simplest + widest)

def prep_preamble():
    sps = int(Fs // CHIP_RATE)
    assert sps == 2, "This script expects 2 samples/chip at 10 MS/s and 5 Mcps."
    base = generate_pn_bpsk(N_CHIPS, sps=sps, seed=7)  # complex64
    # Optional: x = rrc_filter(base, beta=0.25, span=8, sps=sps)
    x = base
    # normalize to <= 0.6 peak amplitude to leave DAC headroom
    x = 0.6 * x / (np.max(np.abs(x)) + 1e-12)
    return x

def make_tx_waveform():
    pre = prep_preamble()
    guard = np.zeros(GUARD_SAMPS, dtype=np.complex64)
    one_burst = np.concatenate([pre, guard])
    # Repeat to build a continuous buffer
    reps = int(Fs*0.1 // len(one_burst)) + 10   # ~100ms worth, repeated
    wf = np.tile(one_burst, reps).astype(np.complex64)
    return wf

def xcorr_fft_upsampled(rx, ref, up=64):
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
    r = np.abs(r)
    r = np.concatenate([r[-(L//2):], r[:L//2]])  # center lags at 0
    return r, L

def quad_peak(y, k):
    if k <= 0 or k >= len(y)-1:
        return float(k)
    y1, y2, y3 = y[k-1], y[k], y[k+1]
    denom = (y1 - 2*y2 + y3)
    delta = 0.5*(y1 - y3)/denom if denom != 0 else 0.0
    return k + delta

def cfarr_threshold(y, win_left=2000, k=6.0):
    # noise stats from a region far from peak; simplistic CFAR
    ref = y[:win_left]
    mu = np.mean(ref); sigma = np.std(ref) + 1e-12
    return mu + k*sigma

def pick_leading_edge(r, up=UPS):
    # find first crossing then local-max refine
    T = cfarr_threshold(r, win_left=min(4000, len(r)//8), k=6.0)
    # search around center (lag 0)
    center = len(r)//2
    search = r[center- int(0.01*len(r)) : center + int(0.01*len(r))]  # ±1% window
    idx0 = np.argmax(search) + (center - int(0.01*len(r)))
    # leading-edge: walk left until below threshold, then pick local max forward
    i = idx0
    while i > 1 and r[i] > T:
        i -= 1
    left = max(i, 1)
    right = min(left + 6*UPS, len(r)-2)  # small window after first-crossing
    kmax = left + int(np.argmax(r[left:right]))
    return quad_peak(r, kmax)

def usrp_make_rx(serial:str):
    u = uhd.usrp.MultiUSRP(f"serial={serial}")
    u.set_clock_source("external")
    u.set_time_source("external")
    u.set_rx_rate(Fs)
    u.set_rx_freq(FC)
    u.set_rx_gain(RX_GAIN)
    st = u.get_rx_stream(uhd.usrp.StreamArgs("fc32","sc16"))
    return USRPHandle(u, st)

def usrp_make_tx(serial:str|None):
    args = "" if serial is None else f"serial={serial}"
    u = uhd.usrp.MultiUSRP(args)
    u.set_clock_source("external")
    u.set_time_source("external")
    u.set_tx_rate(Fs)
    u.set_tx_freq(FC)
    u.set_tx_gain(TX_RF_GAIN)
    st = u.get_tx_stream(uhd.usrp.StreamArgs("fc32","sc16"))
    return USRPHandle(u, st)

def sync_all_to_pps(devs):
    for d in devs:
        d.u.set_time_next_pps(uhd.types.TimeSpec(0.0))
    time.sleep(1.1)  # wait for PPS tick

def start_rx_aligned(devs, start_in_sec=0.5, nsamps:int=0):
    t0 = devs[0].u.get_time_now() + uhd.types.TimeSpec(start_in_sec)
    cmd = uhd.types.StreamCmd(uhd.types.StreamMode.START_CONT)
    cmd.stream_now = False
    cmd.time_spec = t0
    for d in devs:
        d.st.issue_stream_cmd(cmd)
    # allocate buffers
    num = int(nsamps) if nsamps>0 else int(Fs*LOCATE_DURATION_S)
    bufs = [np.zeros(num, dtype=np.complex64) for _ in devs]
    mds = [uhd.types.RXMetadata() for _ in devs]
    # receive once (blocking). For robustness you can loop/accumulate frames.
    for i, d in enumerate(devs):
        n = d.st.recv([bufs[i]], mds[i], timeout=3.0)
        if n != len(bufs[i]):
            print(f"[WARN] RX {i} expected {len(bufs[i])} got {n}")
    # stop
    stop = uhd.types.StreamCmd(uhd.types.StreamMode.STOP_CONT)
    for d in devs:
        d.st.issue_stream_cmd(stop)
    return bufs

def compute_toa_seconds(iq, preamble, up=UPS):
    r, L = xcorr_fft_upsampled(iq, preamble, up=up)
    k_hat = pick_leading_edge(r, up=up)
    # lag index to time: (k - L/2) / (Fs * up)
    tau = (k_hat - L/2) / (Fs * up)
    return tau

def multilateration_2d(anchors: np.ndarray, tdoa_sec: np.ndarray, ref:int=0):
    # Gauss-Newton on ranges difference
    x = np.mean(anchors, axis=0).astype(float)
    for _ in range(20):
        H = []
        r = []
        for i in range(len(anchors)):
            if i == ref: continue
            di = np.linalg.norm(x - anchors[i]) + 1e-9
            dr = np.linalg.norm(x - anchors[ref]) + 1e-9
            model = di - dr
            meas  = c * tdoa_sec[i]
            r.append(model - meas)
            gi = (x - anchors[i])/di - (x - anchors[ref])/dr
            H.append(gi)
        H = np.vstack(H); r = np.array(r)
        dx, *_ = np.linalg.lstsq(H, -r, rcond=None)
        x = x + dx
        if np.linalg.norm(dx) < 1e-3: break
    return x

def run_tx():
    print("[TX] starting…")
    txh = usrp_make_tx(TX_SERIAL)
    sync_all_to_pps([txh])
    wf = make_tx_waveform()
    # Start at next 0.5s
    t0 = txh.u.get_time_now() + uhd.types.TimeSpec(0.5)
    cmd = uhd.types.StreamCmd(uhd.types.StreamMode.START_CONT)
    cmd.stream_now = False
    cmd.time_spec = t0
    # UHD TX: send forever in a loop
    print("[TX] streaming at", t0.get_real_secs(), "seconds")
    time.sleep(0.6)
    md = uhd.types.TXMetadata()
    md.start_of_burst = True
    md.end_of_burst = False
    md.has_time_spec = True
    md.time_spec = t0
    while True:
        sent = txh.st.send([wf], md)
        md.start_of_burst = False
        md.has_time_spec = False
        # throttle a bit (waveform is ~100 ms; UHD buffers usually handle this)
        # break with Ctrl+C

def run_calibration():
    print("[CAL] building receivers…")
    devs = [usrp_make_rx(s) for s in RX_SERIALS]
    sync_all_to_pps(devs)
    nsamps = int(Fs*CALIB_DURATION_S)
    print(f"[CAL] capturing {nsamps} samples on {len(devs)} receivers…")
    bufs = start_rx_aligned(devs, start_in_sec=0.5, nsamps=nsamps)
    pre = prep_preamble()
    toas = [compute_toa_seconds(b, pre, up=UPS) for b in bufs]
    # biases relative to receiver 0 (seconds)
    ref = 0
    biases = [float(toas[i] - toas[ref]) for i in range(len(toas))]
    ns_biases = [b*1e9 for b in biases]
    print("[CAL] per-receiver biases relative to RX0 (ns):", ns_biases)
    data = {"serials": RX_SERIALS, "bias_sec": biases}
    Path(BIAS_FILE).write_text(json.dumps(data, indent=2))
    print(f"[CAL] saved {BIAS_FILE}")

def run_locate():
    print("[LOC] loading biases…")
    if not Path(BIAS_FILE).exists():
        raise SystemExit("Bias file not found. Run 'calibrate' first.")
    data = json.loads(Path(BIAS_FILE).read_text())
    assert data["serials"] == RX_SERIALS, "Receiver order changed; re-run calibration."
    bias = np.array(data["bias_sec"], dtype=float)
    devs = [usrp_make_rx(s) for s in RX_SERIALS]
    sync_all_to_pps(devs)
    nsamps = int(Fs*LOCATE_DURATION_S)
    print(f"[LOC] capturing {nsamps} samples on {len(devs)} receivers…")
    bufs = start_rx_aligned(devs, start_in_sec=0.5, nsamps=nsamps)
    pre = prep_preamble()
    toas = np.array([compute_toa_seconds(b, pre, up=UPS) for b in bufs], dtype=float)
    # apply bias correction (remove per-RX fixed delays)
    toas_corr = toas - bias
    # TDOA relative to reference 0
    tdoa = toas_corr - toas_corr[0]
    print("[LOC] TOAs (ns):", (toas_corr*1e9).round(2).tolist())
    print("[LOC] TDOAs vs RX0 (ns):", (tdoa*1e9).round(2).tolist())
    if len(RX_SERIALS) < 3 or ANCHORS.shape[0] < 3:
        print("[LOC] Need >=3 anchors for 2D fix. Showing hyperbola line only.")
        return
    est = multilateration_2d(ANCHORS, tdoa, ref=0)
    print(f"[LOC] Estimated position (m): x={est[0]:.2f}, y={est[1]:.2f}")
    # quick plot (optional)
    try:
        import matplotlib.pyplot as plt
        xs, ys = ANCHORS[:,0], ANCHORS[:,1]
        plt.figure()
        plt.scatter(xs, ys, marker='^')
        for i,(x,y) in enumerate(ANCHORS):
            plt.text(x+0.2, y+0.2, f"RX{i}")
        plt.scatter([est[0]], [est[1]], marker='x')
        plt.title("TDOA Estimate")
        plt.xlabel("x (m)"); plt.ylabel("y (m)")
        plt.axis('equal'); plt.grid(True)
        plt.savefig("tdoa_estimate.png", dpi=150)
        print("[LOC] Saved tdoa_estimate.png")
    except Exception as e:
        print("[LOC] Plot skipped:", e)

def main():
    ap = argparse.ArgumentParser(description="End-to-end TDOA demo: Tx, Calibrate, Locate")
    ap.add_argument("mode", choices=["tx","calibrate","locate"], help="Run mode")
    args = ap.parse_args()
    if args.mode == "tx":
        run_tx()
    elif args.mode == "calibrate":
        run_calibration()
    elif args.mode == "locate":
        run_locate()

if __name__ == "__main__":
    main()
