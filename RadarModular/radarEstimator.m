% Estimating Radar Parameters
function [outputMat, est] = radarEstimator(data, rP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimating Radar Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rP.D            = data.D;
rP.freq         = data.freq;
TEst            = data.det.TEst;
PWEst           = data.det.PWEst;
% set the waveform (estimated)
[rP] = selectWaveform( PWEst/rP.fs, TEst/rP.fs, 0, true, rP);
[outputMat, est] = radarDetector ( data, rP);

% Calculating  g
[~, ~, g] = majoritySelect (est.estgMF);
est.g=g;

% Calculating error PW
PW = PWEst/rP.fs;   est.PW=PW;

% Calculating  BW
BW = g*PW;  est.BW=BW;

% Calculating  T
T = TEst/rP.fs;  est.T=T;

if rP.debugLog
    fprintf(' Estimated: PW = %.2f (us), PRF = %.2f (KHz),  BW = %.2f (MHz), g = %.2f (GHz/s) \n',  PW/1e-6, 1/T/1e3, BW/1e6, g/1e9 );
end
% Write to results
fprintf(rP.fileID,' Estimated: PW = %.2f (us), PRF = %.2f (KHz),  BW = %.2f (MHz), g = %.2f (GHz/s)  \n', PW/1e-6, 1/T/1e3, BW/1e6, g/1e9 );

[estMode] = selectWaveform( PW, T, BW, false, rP);
est.estMode=estMode;

if rP.debugLog
    fprintf(' Estimated: waveform = %s   \n', estMode );
end
% Write to results
fprintf(rP.fileID,' Estimated: waveform = %s   \n', estMode );

% Estimating Radar CSI
if rP.debugLog
    fprintf(' Estimating Radar RCSI ...  \n' );
end
% Write to results
fprintf(rP.fileID,' Estimating Radar RCSI ...  \n' );

rP.mode = estMode;
%radarCSI( PW, BW, 1/T, data, rP);
if rP.debugLog
    fprintf(' RCSI Saved in JSON file!  \n' );
end
% Write to results
fprintf(rP.fileID,' RCSI Saved in JSON file!  \n' );
%%%%%%%%%%%%%%%%%
end

