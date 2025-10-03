function [cpeaks, TEst, PWEst, curSig, data, rP] = estimateTPW (rxsig, data, fs, Ts, sigScaling, rP, flagSEC)


if ~exist('flagSEC','var')
    flagSEC = false;
end

PLOT_FIG        = false;
DET_DECISION    = false;
NpulseBIN1A     = 400;
NpulseBIN1B     = 600;
NpulseBIN2A     = 4;
NpulseBIN2B     = 4;
NpulseBIN3A     = 2;
NpulseBIN3B     = 3;
rP.percent_peak = 0.4; %0.40; 0.6182;
rP.RadarPresent = true;
rP.numPulses    = 10;
offsetPW        = 0;
totalFreqBins   = 84;
freq            = linspace(0,1, totalFreqBins).';  % Frequency bins
freq            = freq(1:totalFreqBins);
[NN, Nch]       = size(rxsig);
curSig          = rxsig;
eps             = 1e-3;
TT              = 3072 -6; % smallest Inter-pulse interval

[sigBin1A, NpBin1A, PRIBin1A, Ns] = generateSig ('BIN1-A', fs, Ts,  sigScaling);
[sigBin1B, NpBin1B, PRIBin1B, Ns] = generateSig ('BIN1-B', fs, Ts,  sigScaling);
[sigBin2A, NpBin2A, PRIBin2A, Ns] = generateSig ('BIN2-A', fs, Ts,  sigScaling);
[sigBin2B, NpBin2B, PRIBin2B, Ns] = generateSig ('BIN2-B', fs, Ts,  sigScaling);
[sigBin3A, NpBin3A, PRIBin3A, Ns] = generateSig ('BIN3-A', fs, Ts,  sigScaling);
[sigBin3B, NpBin3B, PRIBin3B, Ns] = generateSig ('BIN3-B', fs, Ts,  sigScaling);



%=============   Estimating Inter-pulse Interval TEst ================================
indCh = 1; % selecting channel to perform TEst estimation
rP.percent_peak = 0.4;
vecRx                = rxsig(:,indCh);

vecBin1A = sigBin1A(1:NpBin1A,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin1A, rP.D, rP.freq, false);
    vecBin1A    = sum(abs(outputMat));
end
conVBin1A = abs( conv(vecRx.' , fliplr(vecBin1A') ) ) ;

vecBin1B = sigBin1B(1:NpBin1B,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin1B, rP.D, rP.freq, false);
    vecBin1B    = sum(abs(outputMat));
end
conVBin1B = abs( conv(vecRx.' , fliplr(vecBin1B') ) ) ;

vecBin2A = sigBin2A(1:NpBin2A,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin2A, rP.D, rP.freq, false);
    vecBin2A    = sum(abs(outputMat));
end
conVBin2A = abs( conv(vecRx.' , fliplr(vecBin2A') ) ) ;

vecBin2B = sigBin2B(1:NpBin2B,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin2B, rP.D, rP.freq, false);
    vecBin2B    = sum(abs(outputMat));
end
conVBin2B = abs( conv(vecRx.' , fliplr(vecBin2B') ) ) ;

vecBin3A = sigBin3A(1:NpBin3A,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin3A, rP.D, rP.freq, false);
    vecBin3A    = sum(abs(outputMat));
end
conVBin3A = abs( conv(vecRx.' , fliplr(vecBin3A') ) ) ;

vecBin3B = sigBin3B(1:NpBin3B,1);
if flagSEC
    [outputMat] = SEC_ALG(vecBin3B, rP.D, rP.freq, false);
    vecBin3B    = sum(abs(outputMat));
end
conVBin3B = abs( conv(vecRx.' , fliplr(vecBin3B') ) ) ;

if PLOT_FIG
    figure; plot(real(conVBin1A))
    figure; plot(real(conVBin1B))
    figure; plot(real(conVBin2A))
    figure; plot(real(conVBin2B))
    figure; plot(real(conVBin3A))
    figure; plot(real(conVBin3B))
end

% Making decision which peaks to use given correlation of each waveforms
mB1A = mean(conVBin1A); maxB1A = max(conVBin1A);
mB1B = mean(conVBin1B);  maxB1B = max(conVBin1B);
mB2A = mean(conVBin2A);  maxB2A = max(conVBin2A);
mB2B = mean(conVBin2B);  maxB2B = max(conVBin2B);
mB3A = mean(conVBin3A);  maxB3A = max(conVBin3A);
mB3B = mean(conVBin3B);  maxB3B = max(conVBin3B);

maxAll = [maxB1A, maxB1B maxB2A maxB2B maxB3A maxB3B];
meanAll = [mB1A, mB1B mB2A mB2B mB3A mB3B] ;
indZ = meanAll <eps ;
if sum(indZ) == 6
    meanAll = ones(1, 6 );
elseif sum(indZ) > 0
    minZ = min(meanAll(~indZ));
    meanAll(indZ) = minZ;
end

%[maxV, indV] = max([maxB1A/mB1A, maxB1B/mB1B maxB2A/mB2A maxB2B/mB2B maxB3A/mB3A maxB3B/mB3B]);
[maxV, indV] = max([maxAll./meanAll]);
%----------------------------------------------------------------------

if indV == 1
    rP.numPulses    = NpulseBIN1A+1;
    conVB = conVBin1A;
    NpB  = NpBin1A;
elseif indV == 2
    rP.numPulses    = NpulseBIN1B+1;
    conVB = conVBin1B;
    NpB  = NpBin1B;
elseif indV == 3
    rP.numPulses    = NpulseBIN2A+1;
    conVB = conVBin2A;
     NpB  = NpBin2A;
elseif indV == 4
    rP.numPulses    = NpulseBIN2B+1;
    conVB = conVBin2B;
    NpB  = NpBin2B;
elseif indV == 5
    rP.numPulses    = NpulseBIN3A+1;
    conVB = conVBin3A;
    NpB  = NpBin3A;
elseif indV == 6
    rP.numPulses    = NpulseBIN3B+1;
    conVB = conVBin3B;
    NpB  = NpBin3B;
    
end

[cpeaks, ~, TEst, ~, ~] = findPeaks2 ( conVB.', rP);  % Best estimate of TEst

% Making sure cpeaks smallest inter-pulse interval is larger than the
% smallest T =  3072
[maxPeak, maxPInd] = max(conVB(cpeaks));
cpeaks_ = cpeaks(maxPInd);
sInd = maxPInd;
while (sInd < length(cpeaks))
    sInd = sInd + 1;
    if (cpeaks(sInd) - cpeaks_(end) > TT)
        cpeaks_ = [cpeaks_ cpeaks(sInd)];
    end
end
sInd = maxPInd-1;
while (sInd > 0)
    if (cpeaks_(1) - cpeaks(sInd) > TT)
        cpeaks_ = [cpeaks(sInd) cpeaks_ ];
    end
    sInd = sInd - 1;
end
if cpeaks_>fix(NpB/2)
    cpeaks = cpeaks_-fix(NpB/2); % best estimate of Peaks after MF
else
    cpeaks = cpeaks_;
end

%=============   Estimating Pulse width PWEst ================================

% Having the peaks just select maximum PW and then process it
if 0
    PWEst = 10000;
    [rxsig_, LPW, indPulses] = SliceSignal(cpeaks, PWEst, rxsig, offsetPW);
else
    rxsig_ =  rxsig;
end


TEstAll = []; PWEstAll = []; outputMatAll = [];
for ch = 1: 1%Nch
    [outputRx] = SEC_ALG(rxsig_(:,ch), 16,  freq, false);
    %outputRxAll(:,:, ch) = outputRx;
    %outputRxAll(:, ch) = outputRx(:);
    sumOutPutRx      = sum(abs(outputRx))';
    [rP.percent_peak] = thresholdC(sumOutPutRx);
    [~, ~, TEst_, ~, PWEst_] = findPeaks2 ( sumOutPutRx, rP);
    TEstAll = [TEstAll TEst_]; PWEstAll = [PWEstAll PWEst_];
end

[inxAnomaly, indCh, TEst_] = majoritySelect (TEstAll, 0.5);
[inxA, ~, PWEst]        = majoritySelect (PWEstAll);

%vecRx                = rxsig_(:,indCh);
[~, rP.indCh]        = find(~ismember([1:Nch], inxAnomaly));

mode            = 'OTHER-A';
if  ( PWEst > 35 && PWEst < 90 ) && ( TEst > 2500 && TEst < 3500 )
    mode   = 'BIN1-A';
    NpB    = NpBin1A;
    conVB  = conVBin1A;
    vecB   = vecBin1A;
elseif  ( PWEst > 50 && PWEst < 200 ) && ( TEst > 1500 && TEst < 2500 )
    mode = 'BIN1-B';
    NpB  = NpBin1B;
    conVB = conVBin1B;
    vecB =vecBin1B;
elseif  ( PWEst > 2500 && PWEst < 3600 ) && ( TEst > 250000 && TEst < 350000 )
    mode = 'BIN2-A';
    NpB  = NpBin2A;
    conVB = conVBin2A;
    vecB =vecBin2A;
elseif  ( PWEst > 5500 && PWEst < 6600 ) && ( TEst > 250000 && TEst < 350000 )
    mode = 'BIN2-B';
    NpB  = NpBin2B;
    conVB = conVBin2B;
    vecB =vecBin2B;
elseif  ( PWEst > 450 && PWEst < 750 ) && ( TEst > 500000 && TEst < 700000 )
    mode = 'BIN3-A';
    NpB  = NpBin3A;
    conVB = conVBin3A;
    vecB =vecBin3A;
elseif ( PWEst > 2500 && PWEst < 3700 ) && ( TEst > 350000 && TEst < 450000 )
    mode = 'BIN3-B';
    NpB  = NpBin3B;
    conVB = conVBin3B;
    vecB =vecBin3B;
end % final if


if DET_DECISION && (strcmp(mode, 'OTHER-A')  || ~((TEst < 2*2e-2*rP.fs) && (PWEst < (2*100e-6*rP.fs) )) )
    % we can make a decision that no radar is present since it is way-off
    % from the transmitted radar waveform
    rP.RadarPresent   = false;

    if rP.debugLog
        fprintf(' The radar parameters are way-off !\n');
    end
    % Write to results
    fprintf(rP.fileID,' The radar parameters are way-off !\n');
else
    curSig(:,inxAnomaly) = [];
    [curSig, LPW, indPulses] = SliceSignal(cpeaks, PWEst, curSig, offsetPW);

    % save the pulses
    data.dataDet = curSig;
    if LPW > size(curSig,1)
        LPW = size(curSig,1);
    end
    rP.numPulses = length(cpeaks);
    if rP.numPulses*LPW > size(curSig, 1)
        rP.numPulses = floor(size(curSig, 1)/LPW);
    end

    data.numW        = rP.numPulses;
    data.WINDOW      = LPW;
    data.rWIN        = LPW;

    % save parameters
    data.det.cpeaks = cpeaks;
    data.det.TEst = TEst;
    data.det.PWEst = PWEst;
end

% Estimating Inter-pulse interval
%[cpeaks, peaks, TEst, TEst1, PWEst] = findPeaks2 ( conVB, rP);

end

function [curSig, LPW, indPulses] = SliceSignal(cpeaks, PWEst, curSig, offsetPW)

[NN, ~] = size(curSig);
indPW     = [-offsetPW - fix(PWEst/2):offsetPW + fix(PWEst/2)];
LPW       = length(indPW);
indPulses = (cpeaks)+indPW';
indPulses = indPulses(:);
indPulses = unique(indPulses);
indPulses(indPulses<1) = [];
indPulses(indPulses>NN) = [];
if ~isempty(indPulses)
    curSig = curSig(indPulses, :);
end

end

function [th] = thresholdC(XX)

Xmean = mean(XX);
Xstd = std(XX);
Xmax  = max(XX);
%th    = 1.645536823391266*Xmean+1.914208189073610*Xstd-0.109007794222955*Xmax;
%th    = 1.432*Xmean+0.87*Xstd-0.0198*Xmax;
%th    = 1.1033*Xmean+2.8873*Xstd+0.0752*Xmax-0.8358*Xmax/Xmean;
if Xmean < 1e-3
    th    = 1.3782*Xmean+1.2726*Xstd+0.0431*Xmax+0.0104*Xmax;
else
    th    = 1.3782*Xmean+1.2726*Xstd+0.0431*Xmax+0.00104*Xmax/Xmean;
end
th = th/max(abs(real(XX)));

end


