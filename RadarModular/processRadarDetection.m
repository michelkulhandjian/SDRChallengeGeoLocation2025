% Process radar detection
function  processRadarDetection( NameFile, inputPathDir, outputPathDir, mode, rP)

% ==Radar Detection with different modes
% 1. Simulate - given the waveform, SNR, from the params 
%   a. generateSig()
%   b. transmitRadar()
%   c. addingAWGN()
%   d. Preprocessing() - mode 1.
%   e. radarDetector() 
% 
% 2. IQ dataset - user provides the InputPathDir 
%   a. loadDataset()
%   b. Preprocessing() - mode 1.
%   e. radarDetector()
%   
% 3. SEC dataset - user provides the InputPathDir 
%   a. loadDataset()
%   b. Preprocessing() - mode 3.
%   e. radarDetector()
%
% 4. IQ dataset - user provides the InputPathDir - all IQ samples
%   a. loadDataset()
%   b. Preprocessing() - mode 1.
%   e. radarDetector()

[data, rP] = setParams(NameFile, inputPathDir, outputPathDir, rP);

if rP.debugLog
    fprintf(' Start radar detection process... \n');
end

% Write to results
fprintf(rP.fileID,' Start radar detection process... \n');

% Parameters
fs                = rP.fs;
fc                = rP.fc;
Ts                = rP.Ts; % length of the radar stream in sec
curSNR            = rP.DETSNR;
curWaveF          = rP.DETWaveform;

%% Channel model
chanModel         = rP.chanModel;
isApplyChannel    = true;
sigScaling        = 0.9;  % Scaling of the signal

if mode == 1 || mode == 2
    % Generate the Radar waveform0;
    [sig, Np, PRI, Ns] = generateSig (curWaveF, fs, Ts,  sigScaling);
end

if mode == 1
    % Receiver with Channel enabled by setting it true
    [rxsig] = transmitRadar (sig, chanModel, fs, fc, rP.pPos,  rP.pVel, isApplyChannel);
elseif mode == 2 || mode == 4
    %[rxsig]= processIQtoMAT( [], rP);
    % load IQ dataset and save radar signal in data.sig
    [data, rP] = generate_radar_dataset( rP.inputPathDir, [], 32, rP);
    [rxsig]= data.sig;
    curSig = rxsig;
elseif mode == 3
    [data, rP] = generate_radar_dataset( rP.inputPathDir, [], 64, rP);
end

% Setting parameters
rP.Np        = 60;  % this is for fixed W for processing
if strcmp(curWaveF,'BIN1-A') || strcmp(curWaveF,'BIN1-B')
    rP.numWperPW = 1;
else
    rP.numWperPW = 3;
end
if (mode==1 && isApplyChannel) || (mode==2)
    vec1 = sig(1:Np,1);
    %vec2 = rxsig(1:Np*2,1);
    vec2 = rxsig(1:3000,1);
    conV = abs( conv(vec2.' , fliplr(vec1') ) ) ;
    %[~, ind] = find (conV > max(conV)*0.9);
    [~, ind] = max(conV);
    DelayOffset = ind(1) - Np;
else
    DelayOffset = 0;
end

if (mode == 1 || mode == 2 || mode == 4)
    if (mode == 1 || mode == 2 )
        indSig = [[1:rP.Np*rP.numWperPW]+DelayOffset+[0:PRI:Ns-1]']';
        mSamples = min(numel(indSig), 960);
        rP.indSig = indSig(1:mSamples);
        % Enable this if clear signal testing needed
        %rxsig1 = sig*ones(1,84);
    elseif mode == 4
        rP.indSig = ':';
    end
    if (mode == 1)
        pdBm            = 0; %-113dBm;
        pRMS = rms(rxsig(:))^2;
        rxsig = rxsig/sqrt(pRMS)*10^((pdBm - 30)/20);

        % Add AWGN noise accross all the received antennas mode 1 only
        [curSig, curNoise] = addingAWGN (rxsig, curSNR, 1);
    end
    [data, rP] = preprocessRadar(data, curSig, rP, 1);
end

rP.ispMSs=true;
rP.mode          = 'DETECT';   % LFM Waveform from the Radar Parameters
rP = setRadarRange( rP);
if rP.debugLog
    fprintf(' Starting Radar Detector... \n');
end

% Write to results
fprintf(rP.fileID,' Starting Radar Detector... \n');
[ STFT_s, det] =  radarDetector ( data, rP);

if rP.isSpectrogramF
    rP.Spec       = true;
    [figSTFT_s, ~, ~] = captureSig ( STFT_s, rP, data);

    if ispc
        saveas(figSTFT_s,[rP.imageFolder '\_CombSTFTDet_'  '.png']);
    else
        saveas(figSTFT_s,[rP.imageFolder '/_CombSTFTDet_'  '.png']);
    end
end

save (rP.nameDetRes, 'STFT_s', 'det');

if sum(det.pMSs>rP.DetThreshold)
    if rP.debugLog
        fprintf(' Radar might be present... \n');
    end
    % Write to results
    fprintf(rP.fileID,' Radar might be present... \n');
else
    if rP.debugLog
        fprintf(' No radar is detected... \n');
    end
    % Write to results
    fprintf(rP.fileID,' No radar is detected... \n');
end

fclose(rP.fileID);
end