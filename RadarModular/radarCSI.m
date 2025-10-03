% Estimating Radar CSI
function radarCSI( pwEst, bwEst, prfEst, data, rP)

sampleRate = rP.fs;
pulseSamples = data.WINDOW;
dataMatFilt = data.data.';
capturedPulses = length(dataMatFilt) / pulseSamples;
numAntTot = size(dataMatFilt,1);
%% Segment the the signal up into pulses and take FFT
pulseMat = reshape(dataMatFilt,[numAntTot pulseSamples capturedPulses]);
Nfft = 4096;
fDataMat = zeros(numAntTot,Nfft,capturedPulses);
for pulse = 1:capturedPulses
    indStart = ((pulse -1) * pulseSamples) + 1;
    indEnd = (indStart + pulseSamples) -1;
    fDataMat(:,:,pulse) = fft( pulseMat(:,indStart:indEnd,pulse),Nfft,2);
end

%% distributed detection and parameter estimation
%numAntTot = size(fDataMat,1);
%Nfft = size(fDataMat,2);
numTrials = 10;                                % number of pulses to take from fDataMat
f0 = 0;
Np = fix(pwEst*sampleRate);
[Np, sig] = mySig ( f0, sampleRate, bwEst, pwEst, prfEst, 1, Np);
sig2_f = fft(sig,Nfft);

if numTrials > capturedPulses 
    numTrials = capturedPulses;
end

sig_mat_f = repmat(sig2_f.', numAntTot,1,numTrials);
csi_scs = fDataMat(:,:,1:numTrials)./sig_mat_f;
csi_scs = mean(csi_scs,3);
rcsi_to_json(csi_scs, Nfft, rP.nameRCSIJson, rP.mode, rP.channels);
%write_csi_to_json(csi_scs, Nf, rP.nameRCSIJson, rP.inputPathDir);

end

function [N sig] = mySig ( f0, fs, sweepBW, PW, PRF, numPulses, N)

tau = 0;
Ts = 1/fs;
g1 = sweepBW/PW;
T = fix(1./PRF*fs); % Pulse interval in samples
%Np = fix(PW*fs);
Np = ceil(PW*fs);
nRange = [0:Np-1];
sig = zeros(1, T);
pulse= exp( 1i* ( 2*pi*f0*Ts*(nRange+tau) + pi*g1.*(Ts*(nRange+tau)).^2 ) );
if N > Np
    pulse = [pulse, zeros(1, N-Np)];
    Np = N;
end
sig(1:Np) = pulse;
sig = kron(ones(1,numPulses), sig);
%sig = sig(1:Np).';
sig = sig.';

end



