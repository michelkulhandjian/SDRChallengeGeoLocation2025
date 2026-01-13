
function  crossCorrolateStudies(inputPathDir, outputPathDir, mode, rP)

[data, rP] = setParams([], inputPathDir, outputPathDir, rP);

if rP.debugLog
    fprintf(' Start Cross Correlation Study process... \n');
end

% Write to results
fprintf(rP.fileID,' Start Cross Correlation Study process... \n');

fs              = rP.fs;
rP.crossMode    = 0;       % select the which cross-correlation to use
rP.analyzeCross = true;    % analyse the underlying cross correlation
c               = physconst('LightSpeed');
PLOT_FIG        = false;


if mode == 0
    % load IQ dataset and save signal in data.sig
    [data, rP] = generate_dataset( rP.inputPathDir, [], 32, rP);
else
    % Generate signal based on params.m
    [data, rP] = generate_dataset( rP.inputPathDir, [], 128, rP);
end

[sig]= data.sig;
% Estimate TDOAs via cross-correlation
[tdoaData] = calcDelays(sig, rP);
rangeDifEst=c*tdoaData;

% Compute true TDOA using transmitter and receiving nodes
[tDelay]=calc_Delay(rP.bPos, rP.pPos);
tdoaData = tDelay(2:end)-tDelay(1);
rangeDif =c*tdoaData;

Error = rangeDifEst - rangeDif;
ErrorRMSE = rmse(rangeDifEst,rangeDif, 1);

fprintf("TDOA errors between estimated and true: "); fprintf(" %.2f  ", Error); fprintf('\n');
fprintf("TDOA RMSE   between estimated and true: "); fprintf(" %.2f  ", ErrorRMSE); fprintf('\n');

save (rP.nameGeolocateRes, 'data', 'rP');
fclose(rP.fileID);

end