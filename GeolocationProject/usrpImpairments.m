function [sig] = usrpImpairments(sig, fc, fs, rP, stage, rxIndex)
% usrpImpairments - Apply realistic USRP hardware impairments
%
%   sig = usrpImpairments(sig, fc, fs, rP, 'tx')            B210 TX impairments
%   sig = usrpImpairments(sig, fc, fs, rP, 'rx')            X310 RX impairments
%   sig = usrpImpairments(sig, fc, fs, rP, 'rx_async', idx)  Async RX offsets for receiver idx
%
% Models the following impairments based on AD9361 (B210) and UBX-160 (X310) datasheets:
%   TX:       DAC quantization, IQ imbalance, DC offset, phase noise, PA nonlinearity
%   RX:       IQ imbalance, DC offset, phase noise, ADC quantization
%   RX_ASYNC: Sampling clock offset, carrier freq offset, time offset, phase offset

if ~isfield(rP, 'usrp') || ~rP.usrp.enabled
    return;
end

if ~exist('rxIndex','var')
    rxIndex = 1;
end

cfg = rP.usrp;

switch lower(stage)
    case 'tx'
        sig = applyTxImpairments(sig, fc, fs, cfg);
    case 'rx'
        sig = applyRxImpairments(sig, fc, fs, cfg);
    case 'rx_async'
        sig = applyAsyncRx(sig, fc, fs, cfg, rxIndex);
    otherwise
        error('Stage must be ''tx'', ''rx'', or ''rx_async''');
end

end

%% ========================================================================
%  B210 TX Impairments (AD9361 RF Transceiver)
%  ========================================================================
function sig = applyTxImpairments(sig, fc, fs, cfg)

    % --- 1. DAC Quantization (12-bit) ---
    dacBits = cfg.tx.dacBits;  % 12
    sigMax = max(abs([real(sig(:)); imag(sig(:))]));
    if sigMax > 0
        % Scale to full-scale, quantize, scale back
        scaleFactor = (2^(dacBits-1) - 1) / sigMax;
        I = round(real(sig) * scaleFactor) / scaleFactor;
        Q = round(imag(sig) * scaleFactor) / scaleFactor;
        sig = complex(I, Q);
    end

    % --- 2. TX IQ Imbalance ---
    % AD9361 after TX quad calibration: image ~-50 dB below carrier
    % This corresponds to ~0.3% amplitude imbalance, ~0.2 deg phase imbalance
    gainImbalance_dB = cfg.tx.iqGainImbalance_dB;   % e.g., 0.05 dB
    phaseImbalance_deg = cfg.tx.iqPhaseImbalance_deg; % e.g., 0.3 degrees
    sig = applyIQImbalance(sig, gainImbalance_dB, phaseImbalance_deg);

    % --- 3. TX DC Offset / LO Leakage ---
    % AD9361 spec: -50 dBc at 0 dB attenuation after calibration
    dcOffset_dBc = cfg.tx.dcOffset_dBc;  % e.g., -50
    sigRMS = rms(sig(:));
    if sigRMS > 0
        dcLevel = sigRMS * 10^(dcOffset_dBc/20);
        % Random phase for LO leakage
        rng_state = rng; rng(42, 'twister'); % deterministic for reproducibility
        dcPhase = 2*pi*rand();
        rng(rng_state);
        sig = sig + dcLevel * exp(1i*dcPhase);
    end

    % --- 4. TX Phase Noise ---
    % AD9361 integrated phase noise: ~0.37 deg RMS at 2.4 GHz, ~0.59 deg at 5.5 GHz
    % Interpolate for operating frequency
    phaseNoiseRMS_deg = cfg.tx.phaseNoiseRMS_deg;
    sig = applyPhaseNoise(sig, phaseNoiseRMS_deg, fs);

    % --- 5. TX PA Nonlinearity (3rd-order model) ---
    % AD9361 OIP3: +23 dBm at 800 MHz, +19 dBm at 2.4 GHz, +17 dBm at 5.5 GHz
    if cfg.tx.enablePA
        oip3_dBm = cfg.tx.oip3_dBm;
        txPower_dBm = cfg.tx.outputPower_dBm;
        sig = applyPANonlinearity(sig, oip3_dBm, txPower_dBm);
    end

end

%% ========================================================================
%  X310 RX Impairments (UBX-160 Daughterboard + ADS62P48 ADC)
%  ========================================================================
function sig = applyRxImpairments(sig, fc, fs, cfg)

    % --- 1. RX IQ Imbalance ---
    % UBX-160 / AD9361 RX: ~0.2% gain error, ~0.2 deg phase error
    gainImbalance_dB = cfg.rx.iqGainImbalance_dB;
    phaseImbalance_deg = cfg.rx.iqPhaseImbalance_deg;
    sig = applyIQImbalance(sig, gainImbalance_dB, phaseImbalance_deg);

    % --- 2. RX DC Offset ---
    % Residual DC after calibration
    dcOffset_dBc = cfg.rx.dcOffset_dBc;
    sigRMS = rms(sig(:));
    if sigRMS > 0
        dcLevel = sigRMS * 10^(dcOffset_dBc/20);
        rng_state = rng; rng(137, 'twister');
        dcPhase = 2*pi*rand();
        rng(rng_state);
        sig = sig + dcLevel * exp(1i*dcPhase);
    end

    % --- 3. RX Phase Noise ---
    phaseNoiseRMS_deg = cfg.rx.phaseNoiseRMS_deg;
    sig = applyPhaseNoise(sig, phaseNoiseRMS_deg, fs);

    % --- 4. ADC Quantization (14-bit, ADS62P48) ---
    adcBits = cfg.rx.adcBits;  % 14
    sigMax = max(abs([real(sig(:)); imag(sig(:))]));
    if sigMax > 0
        scaleFactor = (2^(adcBits-1) - 1) / sigMax;
        I = round(real(sig) * scaleFactor) / scaleFactor;
        Q = round(imag(sig) * scaleFactor) / scaleFactor;
        sig = complex(I, Q);
    end

end

%% ========================================================================
%  Asynchronous Receiver Offsets (independent clocks per X310)
%  ========================================================================
function sig = applyAsyncRx(sig, fc, fs, cfg, rxIndex)
% Apply per-receiver asynchronous offsets when each X310 has its own clock.
% Uses rxIndex as a deterministic seed so each receiver gets a unique but
% reproducible set of offsets.

    if ~cfg.rx.asyncEnabled
        return;
    end

    N = length(sig);
    t = (0:N-1).' / fs;  % time vector (seconds)

    % Use rxIndex-based seed for reproducible per-receiver randomness
    rng_state = rng;
    rng(1000 + rxIndex, 'twister');

    % --- 1. Sampling Clock Offset (SCO) ---
    % Each X310 TCXO has up to +/-2.5 ppm offset from nominal.
    % This causes a time-varying resampling: effective sample rate = fs*(1 + ppm*1e-6)
    % Modeled as a linear phase ramp (fractional delay that grows over time).
    clockOffset_ppm = cfg.rx.clockOffset_ppm * (2*rand() - 1);  % uniform in [-max, +max]
    if clockOffset_ppm ~= 0
        % Fractional sample drift per sample = ppm * 1e-6
        % Total drift over N samples = N * ppm * 1e-6 samples
        % Apply via linear interpolation (resample)
        driftPerSample = clockOffset_ppm * 1e-6;
        % New sample indices (slightly stretched/compressed time axis)
        newIdx = (0:N-1).' * (1 + driftPerSample);
        % Interpolate I and Q separately using cubic interpolation
        origIdx = (0:N-1).';
        I_resampled = interp1(origIdx, real(sig), newIdx, 'pchip', 0);
        Q_resampled = interp1(origIdx, imag(sig), newIdx, 'pchip', 0);
        sig = complex(I_resampled, Q_resampled);
    end

    % --- 2. Carrier Frequency Offset (CFO) ---
    % Each X310 LO has a small frequency error relative to the TX LO.
    % This manifests as a linearly rotating phase: exp(j*2*pi*df*t)
    freqOffset_Hz = cfg.rx.freqOffset_Hz * (2*rand() - 1);  % uniform in [-max, +max]
    if freqOffset_Hz ~= 0
        sig = sig .* exp(1i * 2 * pi * freqOffset_Hz * t);
    end

    % --- 3. Random Time Offset (sample-level) ---
    % Each receiver starts capturing at a slightly different time.
    % Modeled as a circular shift of the signal.
    maxDelay = cfg.rx.maxTimeOffset_samples;
    timeOffset = randi([0, maxDelay]);
    if timeOffset > 0
        sig = circshift(sig, timeOffset);
        % Zero out the wrapped portion (causal — signal hasn't arrived yet)
        sig(1:timeOffset) = 0;
    end

    % --- 4. Random Initial Phase Offset ---
    % Each X310 LO starts at a random phase since there is no shared PPS.
    if cfg.rx.randomPhaseOffset
        phiOffset = 2*pi*rand();
        sig = sig * exp(1i * phiOffset);
    end

    % Restore RNG state
    rng(rng_state);

end

%% ========================================================================
%  Helper Functions
%  ========================================================================

function sig = applyIQImbalance(sig, gainImbalance_dB, phaseImbalance_deg)
% Apply IQ gain and phase imbalance
%   gainImbalance_dB: amplitude imbalance between I and Q (dB)
%   phaseImbalance_deg: phase error between I and Q (degrees)

    g = 10^(gainImbalance_dB/20);  % linear gain imbalance
    phi = phaseImbalance_deg * pi/180;  % phase imbalance in radians

    % IQ imbalance matrix model:
    %   I' = I
    %   Q' = g * (sin(phi)*I + cos(phi)*Q)
    I = real(sig);
    Q = imag(sig);
    Q_impaired = g * (sin(phi)*I + cos(phi)*Q);
    sig = complex(I, Q_impaired);
end

function sig = applyPhaseNoise(sig, phaseNoiseRMS_deg, fs)
% Apply phase noise to signal
%   phaseNoiseRMS_deg: RMS phase noise in degrees
%   fs: sample rate (used for noise bandwidth shaping)

    if phaseNoiseRMS_deg == 0
        return;
    end

    N = length(sig);
    phaseNoiseRMS_rad = phaseNoiseRMS_deg * pi/180;

    % Generate 1/f-shaped phase noise (more realistic than white phase noise)
    % PLL phase noise has ~1/f^2 characteristic close-in, flattening at higher offsets
    f = (0:N-1)'/N * fs;
    f(1) = f(2);  % avoid division by zero

    % Simple Lorentzian-like PSD: flat below corner, 1/f^2 above
    f_corner = 1e3;  % PLL loop bandwidth ~1 kHz
    H = 1 ./ sqrt(1 + (f/f_corner).^2);
    H = H / rms(H);  % normalize

    % Generate white noise, shape it, scale to desired RMS
    whiteNoise = randn(N, 1);
    shapedNoise = real(ifft(fft(whiteNoise) .* H));
    shapedNoise = shapedNoise / rms(shapedNoise) * phaseNoiseRMS_rad;

    % Apply as phase modulation
    sig = sig .* exp(1i * shapedNoise);
end

function sig = applyPANonlinearity(sig, oip3_dBm, txPower_dBm)
% Apply 3rd-order nonlinearity model (memoryless Saleh-like)
%   oip3_dBm: output 3rd-order intercept point in dBm
%   txPower_dBm: average TX output power in dBm
%
% The AM/AM model:  y = a1*x + a3*x*|x|^2
% where a3 is derived from OIP3

    % Convert to linear power
    P_oip3 = 10^(oip3_dBm/10) * 1e-3;  % Watts
    P_tx = 10^(txPower_dBm/10) * 1e-3;  % Watts

    % Normalize signal to match TX power level
    sigPower = mean(abs(sig(:)).^2);
    if sigPower == 0
        return;
    end
    sig = sig * sqrt(P_tx / sigPower);

    % For a 3rd-order model: OIP3 = sqrt(4*a1^3 / (3*|a3|))
    % With a1 = 1 (linear gain normalized out): |a3| = 4/(3*P_oip3)
    a1 = 1;
    a3 = -4 / (3 * P_oip3);  % negative for compression

    % Apply AM/AM
    sig = a1 * sig + a3 * sig .* abs(sig).^2;

    % Restore original power scaling
    newPower = mean(abs(sig(:)).^2);
    if newPower > 0
        sig = sig * sqrt(sigPower / newPower);
    end
end
