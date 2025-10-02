function [sig, Np, PRI, Ns] = generateSig (mode, fs, Ts, bw, sigScaling, rP)


switch mode
    case "PCW"
        % Configure phase-coded waveform
        N = 1024;                                        % Number of subbands
        numChip = 2^nextpow2(N)-1;                       % Number of chips in one phase-coded waveform
        tchip = 1/bw;                                    % Chip duration (s)
        tWave = numChip * tchip;                         % Modulation period (s)
        prf = 1/tWave;                                   % PRF (1/s)
        Np = 0; PRI = 0; Ns = 0;

        % Configure the phase coded waveform as the maximum length sequence
        pmcwWaveform = phased.PhaseCodedWaveform('Code','Maximum Length Sequence', ...
            'SampleRate',fs,'NumChips',numChip,'ChipWidth',tchip,'PRF',prf);
        sig = pmcwWaveform();
    %case "QPSK"
    case {'BPSK', 'QPSK', '8PSK', 'D8PSK', '16PSK', 'QAM16', 'QAM32', 'QAM64', 'QAM256', 'NXDN48', 'NXDN96', 'TONE', 'PAM2', 'PAM4', 'PAM8', 'PAM16'} 
        Np = 0; PRI = 0; Ns = 0;
        [IQout,  bufLen, fc, packet_len] = gen_signal ( mode, bw, fs, Ts, rP); 
        sig = IQout.';
    otherwise
        disp("Option Radar");
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


end



