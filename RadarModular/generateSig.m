function [sig, Np, PRI, Ns] = generateSig (mode, fs, Ts,  sigScaling)
    [rP] = radarParams(mode);
    numPulses = ceil(Ts*rP.prf);
    Ns        = ceil(Ts*fs);
    Np        = ceil(rP.pulseWidth*fs);
    PRI       = fix(1./rP.prf*fs);
    %lfmSigGen  = phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',rP.pulseWidth,...
    %    'SweepBandwidth',rP.sweepBW,'PRF',rP.prf,'NumPulses',numPulses);
    %sig        = lfmSigGen(); % Generated radar Signal
    [~, sig] = mySig ( 0, fs, rP.sweepBW, rP.pulseWidth, rP.prf, numPulses, Np );
    sig        = sig(1:Ns)*sigScaling;
end



