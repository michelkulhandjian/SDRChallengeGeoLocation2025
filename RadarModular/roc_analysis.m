% Radar Detection Performance using ROC analysis
function roc_analysis( inputPathDir, outputPathDir, mode, rP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROC analysis for different number of SNR and Ant integration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: for ROC, we can work on IQ samples only, and filtering is required.
% We have two modes of operations
% 1. Simulate - given the waveform, SNRs, reporting window (W), P number of antennas per sub-band, sparsness from the params 
%   a. generateSig()
%   b. transmitRadar()
%   c. addingAWGN()
%   d. Preprocessing() - mode 1.
%   e. RadarDetection()
%   f. generateROC() 
%   
% 2. IQ dataset - given the waveform, SNRs, reporting window (W), P number of antennas per sub-band, sparsness from the params 
%   a. loadDataset()
%   b. addingAWGN()
%   c. Preprocessing() - mode 1.
%   d. RadarDetection()
%   e. generateROC() 

[data, rP] = setParams([], inputPathDir, outputPathDir, rP);

if rP.debugLog
    fprintf(' Start ROC Analysis Process for %s waveform... \n', rP.ROCWaveform);
end
% Write to results
fprintf(rP.fileID,' Start ROC Analysis Process for %s waveform... \n', rP.ROCWaveform);

% Results to be save in this file
nameDetRes        = rP.nameDetRes;

% Scenarios to plot/present
% Parameters we are tunning for
SNR               = rP.ROCsnr;
ROCnAntsperGrp    = rP.ROCnAntsperGrp;
ROCWindow         = rP.ROCWindow;
ROCNspelem        = rP.ROCNspelem;
RWIN              = ROCWindow;
curWaveF          = rP.ROCWaveform;
ROCFig            = rP.ROCFig;
Nframe            = 6;

% Parameters
fs                = rP.fs;
fc                = rP.fc;
Ts                = rP.Ts; % length of the radar stream in sec
WAVEFORMS         = rP.WAVEFORMS;

%% Channel model
chanModel         = rP.chanModel;
isApplyChannel    = true;
sigScaling        = 0.9;  % Scaling of the signal

% save it to res structure
res(1).SNR            = SNR;
res(1).ROCnAntsperGrp = ROCnAntsperGrp;
res(1).ROCWindow      = ROCWindow;
res(1).ROCNspelem     = ROCNspelem;
res(1).ROCWaveform    = curWaveF;
res(1).ROCFig         = ROCFig;

% Generate the Radar waveform0;
[sig, Np, PRI, Ns] = generateSig (curWaveF, fs, Ts,  sigScaling);

if mode == 1
    % Receiver with Channel enabled by setting it true
    [rxsig] = transmitRadar (sig, chanModel, fs, fc, rP.pPos,  rP.pVel, isApplyChannel);
elseif mode == 2
    %[rxsig]= processIQtoMAT( [], rP);
    % load IQ dataset and save radar signal in data.sig
    [data, rP] = generate_radar_dataset( rP.inputPathDir, [], 32, rP);
    [rxsig]= data.sig;
end

% Setting parameters
rP.Np        = 60;  % this is for fixed W for processing
if strcmp(curWaveF,'BIN1-A') || strcmp(curWaveF,'BIN1-B')
    rP.numWperPW = 1;
else
    rP.numWperPW = 3;
end
if isApplyChannel
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
indSig = [[1:rP.Np*rP.numWperPW]+DelayOffset+[0:PRI:Ns-1]']';
mSamples = min(numel(indSig), 960);
rP.indSig = indSig(1:mSamples); 
% Enable this if clear signal testing needed
%rxsig1 = sig*ones(1,84);

pdBm            = 0; %-113dBm;
pRMS = rms(rxsig(:))^2;
%pRMS = rms(rxsig(rP.indSig,1))^2;
%pRMS = 100;
rxsig = rxsig/sqrt(pRMS)*10^((pdBm - 30)/20);

Nsnr              = length(SNR);
NantI             = length(ROCnAntsperGrp);
Nwin              = length(ROCWindow);
Nsp               = length(ROCNspelem);
rP.ispMSs=true;
rP.mode          = 'DETECT';   % LFM Waveform from the Radar Parameters
rP = setRadarRange( rP);
rP.label  = [];
% Start Sweeping SNR
for indSNR = 1:Nsnr
    curSNR = SNR(indSNR);
    rP.label{indSNR} = int2str(curSNR);
    % Start Sweeping AntIntperGrp
    for indAntInt = 1:NantI
        %  Antenna Integration
        data.nAntsperGrp  = ROCnAntsperGrp(indAntInt);
        % Start Sweeping reporting window
        for indW = 1:Nwin
            data.WINDOW = ROCWindow(indW);
            % Start Sweeping sparseness
            for indSP = 1: Nsp
                data.num_sparse_elem = ROCNspelem(indSP);
                %====
                h1MSs = []; h0MSs = [];

                for indFrame = 1:Nframe
                    fprintf('SNR: %d/%d, AntInt: %d/%d, R.Window: %d/%d, SP: %d/%d, FRAME: %d/%d \n', ...
                        indSNR, Nsnr, ...
                        indAntInt, NantI, ...
                        indW, Nwin, ...
                        indSP, Nsp, ...
                        indFrame, Nframe);
                    if false
                        curSig            = awgn(rxsig, curSNR,'measured');
                        curNoise          = curSig - rxsig;
                    else
                        % Add AWGN noise accross all the received antennas
                        [curSig, curNoise] = addingAWGN (rxsig, curSNR, mode);
                    end
                    %curSig = rxsig1;
                    [data, rP] = preprocessRadar(data, curSig, rP, 1);
                    [ STFT_s, detS] =  radarDetector ( data, rP);
                    [data, rP] = preprocessRadar(data, curNoise, rP, 1);
                    [ STFT_sN, detN] =  radarDetector ( data, rP);

                    h1MSs=[h1MSs detS.pMSs]; h0MSs=[h0MSs detN.pMSs];
                end
                %%%
                if rP.ROCFig  == 7 || rP.ROCFig  == 8
                    res(indSNR).sigpMSs{indAntInt}{indSP} = h1MSs; %SNR, (AntPerGrp,W)
                    res(indSNR).noisepMSs{indAntInt}{indSP} = h0MSs;
                    [res(indSNR).pd{indAntInt}{indSP}, res(indSNR).pfa{indAntInt}{indSP}] = roc_curve(h1MSs, h0MSs);
                else

                    res(indSNR).sigpMSs{indAntInt}{indW} = h1MSs; %SNR, (AntPerGrp,W)
                    res(indSNR).noisepMSs{indAntInt}{indW} = h0MSs;
                    [res(indSNR).pd{indAntInt}{indW}, res(indSNR).pfa{indAntInt}{indW}] = roc_curve(h1MSs, h0MSs);
                end

                save (nameDetRes, 'res', 'rP');

            end

        end
    end
end

%plotROC ( res, rP )

fclose(rP.fileID);
end

% ROC Analysis
function [sim_pd, sim_pfa] = roc_curve ( h1, h0)

if isempty(h1) | isempty(h0)
    sim_pd = 0;
    sim_pfa = 0;
    return;
end
select_histogram_plot = 0;
select_roc_plot = 0;

L = length(h1);
h1a = abs(h1);
h0a = abs(h0);

thresh_low = min([h1a,h0a]);
thresh_hi  = max([h1a,h0a]);
nbins = 100;

binedges = linspace(thresh_low,thresh_hi,nbins);
if select_histogram_plot == 1
    figure
    histogram(h0a,binedges)
    hold on
    histogram(h1a,binedges)
    hold off
    title('Target-Absent Vs Target-Present Histograms')
    legend('Target Absent','Target Present')
end

nbins = 1000;
%fprintf( 'thresh_low %d thresh_hi %d sec \n', thresh_low, thresh_hi);
thresh_steps = linspace(thresh_low,thresh_hi,nbins);
sim_pd = zeros(1,nbins);
sim_pfa = zeros(1,nbins);
for k = 1:nbins
    thresh = thresh_steps(k);
    sim_pd(k) = sum(h1a >= thresh);
    sim_pfa(k) = sum(h0a >= thresh);
end
sim_pd = sim_pd/L;
sim_pfa = sim_pfa/L;


pfa_diff = diff(sim_pfa);
idx = (pfa_diff == 0);
sim_pfa(idx) = [];
sim_pd(idx) = [];

minpfa = 1e-6;
N = sum(sim_pfa >= minpfa);
sim_pfa = fliplr(sim_pfa(1:N)).';
sim_pd = fliplr(sim_pd(1:N)).';

if select_roc_plot  == 1
    %semilogx(sim_pfa,sim_pd,'r.')
    figure
    semilogx(sim_pfa,sim_pd,'r')
    hold on
    title(' ROC Curves')
    xlabel('Pfa')
    ylabel('Pd')
    grid on
    legend('Proposed','Location','SE')
end
end
