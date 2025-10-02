function [Pxx,F] = plot_psd(x, fs, plotFlag, varargin)
 
% This function estimates the spectral density of the input signal.
%
% [Pxx,F] = plot_psd(x, fs, plotFlag, varargin)
%
% OUTPUT parameters:
% Pxx NormalisedPower spectral density of the signal
% F frequency vector
%
% INPUT parameters:
% x Input signal
% fs Sampling frequency of input signal (optional, default = 1).
% plotFlag Plot the PSD. This argument is optional. The default is to plot if
% no output arguments are used, and not to plot if output arguments
% are used.
% Any more input arguments are passed to plot.

 
 
if ~exist('fs', 'var')
fs = 1;
end
 
L = length(x);
% set up parameter for PSD
NFFT = 2048/2;
if L < NFFT
    NFFT = 2^(nextpow2(L)-1);
end
%WINDOW = chebwin(NFFT,100); % Filter BW = 1.94 * fs/NFFT
%WINDOW = hanning(NFFT); % Filter BW = 1.5 * fs/NFFT
WINDOW = kaiser(NFFT,10); % Filter BW = 1.85 * fs/NFFT
NOVERLAP = NFFT/2;
 
% Calculate PSD
[Pxx,F] = pwelch(x, WINDOW, NOVERLAP, NFFT, fs, 'twosided');

% Shift frequencies
Pxx = abs(fftshift(Pxx'));
ind = find(F>= fs/2);
F(ind) = F(ind) - fs;
F = fftshift(F);
 
% Plot if wanted
if ~exist('plotFlag', 'var')
if nargout == 0
plotFlag = 1;
else
plotFlag = 0;
end
end
 
if plotFlag
plot(F,10*log10(Pxx), varargin{:});
end
