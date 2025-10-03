% Process Radar signals for estimation
function  processRadarEstimation( NameFile, inputPathDir, outputPathDir, mode, rP)

[data, rP] = setParams(NameFile, inputPathDir, outputPathDir, rP);

if rP.debugLog
    fprintf(' Start radar estimation process... \n');
end

% Write to results
fprintf(rP.fileID,' Start radar estimation process... \n');

% Parameters
fs                = rP.fs;
fc                = rP.fc;
Ts                = rP.Ts; % length of the radar stream in sec
curSNR            = rP.DETSNR;
curWaveF          = rP.DETWaveform;
rP.RadarPresent   = true;

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
    [data, rP] = generate_radar_dataset( rP.inputPathDir, rP.outputPathDir, 32, rP);
    [rxsig]= data.sig;
    curSig = rxsig;
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
        % Estimating T and PW
        [cpeaks, TEst, PWEst, curSig, data, rP] = estimateTPW (rxsig, data, fs, Ts, sigScaling, rP);
    end

    if (mode == 1)
        pdBm            = 0; %-113dBm;
        pRMS = rms(rxsig(:))^2;
        rxsig = rxsig/sqrt(pRMS)*10^((pdBm - 30)/20);

        % Add AWGN noise accross all the received antennas mode 1 only
        [curSig, curNoise] = addingAWGN (rxsig, curSNR, 1);
    end
    if rP.RadarPresent
        [data, rP] = preprocessRadar(data, curSig, rP, 4);
    end
end

if rP.RadarPresent
    rP.ispMS=true; rP.mode = 'ESTIMATE';
    rP = setRadarRange( rP);
    if rP.debugLog
        fprintf(' Starting Radar Parameter estimation... \n');
    end
    % Write to results
    fprintf(rP.fileID,' Starting Radar Parameter estimation... \n');
    % Here we can pass the course estimated parameters to our fine
    % estimator.
    %data.det = det;
    rP.detResponse = true; % at this point we are not passing it to Estimatoor
    [outputMatEst, est] = radarEstimator(data, rP);

    if rP.isSpectrogramF
        rP.Spec       = true;
        [figoutputEst, ~, ~] = captureSig ( outputMatEst, rP, data);

        if ispc
            saveas(figoutputEst,[rP.imageFolder '\_CombSECEst_'  '.png']);
        else
            saveas(figoutputEst,[rP.imageFolder '/_CombSECEst_'  '.png']);
        end
    end
    save (rP.nameEstRes, 'outputMatEst',  'est');
end

if rP.debugLog
    fprintf(' Completed Estimation!  \n' );
end
% Write to results
fprintf(rP.fileID,' Completed Estimation!  \n' );

fclose(rP.fileID);
end