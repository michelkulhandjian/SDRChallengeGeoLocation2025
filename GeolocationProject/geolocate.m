
function  geolocate(NameFile, inputPathDir, outputPathDir, mode, rP)

% Geolocate the transmitter position for the given information and parameters set at params.m

% mode composed of bits: |_7_|_6_|_5_|_4_|_3_|_2_|_1_|_0_|
% 0-bit (1)   if set load from IQ for the given parameter settings and correlate to obtain TDOAs
% 1-bit (2)   if set generate IQ and use correlation to obtain TDOAs
% 2-bit (4)   if set compute true TDOAs 
% 3-bit (8)   if set read TDOAs from tdoa.txt
% 4-bit (16)  if set read from cvs file provided in the input folder
% 5-bit (32)  if set use use LAT/LONG if not set then use Euclidean coordinates 
% 6-bit (64)  if set provide plots for different SNR provided at params.m if not set geolocate the position


[data, rP] = setParams(NameFile, inputPathDir, outputPathDir, rP);

if rP.debugLog
    fprintf(' Start geolocation process... \n');
end

% Write to results
fprintf(rP.fileID,' Start geolocation process... \n');

rP.isEuclidean  = ~bitand(mode , 32);
rP.isALg1       = false;
rP.isNonLinear  = false;
fs              = rP.fs;
rP.crossMode    = 0;       % select the which cross-correlation to use
rP.analyzeCross = true;    % analyse the underlying cross correlation
c               = physconst('LightSpeed');
PLOT_FIG        = false;

%% Alg1
if rP.isALg1
    %rP.SampleRateSource  = true;
    rP.Measurement = 'TDOA';
    rP.SampleRate  = fs;
    rP.NumEstimates = 1;
    rP.DelayOffsetInputPort = false;
    rP.VarianceOutputPort  = true;
    rP.SampleRateSource = 'Input port';
    rP.NoisePowerSource = 'Property';
    rP.NoisePower  = 1.0;
    rP.FFTCross    = true;
    rP.InterPolCross = false;
end


if bitand(mode , 3)
    if bitand(mode , 1)
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

    %% Alg1
    if rP.isALg1 
        [Y, rP, var] = tdoaEst(rP, sig,fs);
        rP.Y = Y; rP.var = var;

        if PLOT_FIG
            figure
            plotTDOASpectrumM(rP,'AnchorPairIndex',1,'MaxDelay',1000e-9);
            %plotTOASpectrumM(rP, freqSpacing,'AnchorIndex',1,'MaxDelay',500e-9);
        end
    end
elseif bitand(mode , 4)
    % Compute true TDOA using transmitter and receiving nodes
    [tDelay]=calc_Delay(rP.bPos, rP.pPos);
    tdoaData = tDelay(2:end)-tDelay(1);
    rangeDif =c*tdoaData;
    rP.Y = tdoaData; rP.var = 0.00002*abs(randn(size(tdoaData)));
elseif bitand(mode , 8)
    fileTdoa = 'tdoa4.txt';
    %opts = "'decimalSeparator',','";
    data  = readmatrix(fileTdoa);

    if rP.isEuclidean
        tdoaData = data;
    else
        toa = data(:,3)'/c;
        geoData  = data(:, 1:end-1)';
        rP.geoData  = geoData;
        [baseStation] = geoToEnu (geoData, 2);
        rP.bPos = baseStation;
        tdoaData = toa(2:end)-toa(1);
        
    end
    rP.Y = tdoaData; rP.var = 0.0001*sqrt(mean(tdoaData.^2))*abs(randn(size(tdoaData)));
elseif bitand(mode , 16)
    fileTdoa = 'tdoa2.csv';
    tdoaData = readmatrix(fileTdoa); %
    rP.mode = 'GeolocateMultiplePosition';
end

rP.tdoaData = tdoaData;
if length(tdoaData)==2
    rP.tdoa32 = tdoaData(2)- tdoaData(1);
end

if bitand(mode , 64)
    rP.mode = 'GeolocateSimulation';
else
    rP.mode = 'GeolocatePosition';
end


if rP.debugLog
    fprintf(' Starting Localizer... \n');
end

% Write to results
fprintf(rP.fileID,' Starting Localizer... \n');

[rP] = localize2 ( data, rP);

save (rP.nameGeolocateRes, 'data', 'rP');
fclose(rP.fileID);

end