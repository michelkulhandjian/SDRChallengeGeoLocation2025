function runProcesses (inputPathDir, outputPathDir , runType, mode )

%close all; clear all;

% Enable debug logs
rP.debugLog       = true;

if ~exist('runType','var')
    runType   = "ESTIMATION";    % Options: {"GENDATASET", "DETECTION", "ROC", "ESTIMATION", "ESTIMATION_ANALYSIS"}
end

if ~exist('mode','var')
    mode = 4;  % 3 for generating radar and saving the dataset (.dat)/ 4 is for Estimation given the Dataset .dat format
end

if rP.debugLog
    fprintf( 'Running as runType: %s with mode %d \n', runType, mode);
end

if ~exist('inputPathDir','var') ||  (exist('inputPathDir', 'var') && isempty(inputPathDir))
    %inputPathDir  = 'C:\Users\michel.kulhandjian\Desktop\Projects\RadarProject\Code\Matlab\radarDetection\DrondeData\data\Oscar\BIN1-A';
    %inputPathDir  = 'C:\Users\michel.kulhandjian\Desktop\Projects\RadarProject\Code\Matlab\radarDetection\DrondeData\data\Oscar\gnb\noise';
    inputPathDir  = pwd;
     %inputPathDir               = 'C:\Users\michel.kulhandjian\Desktop\ground1';
end
if ~exist('outputPathDir','var') || (exist('outputPathDir', 'var') && isempty(outputPathDir))
    outputPathDir  = 'C:\Users\michel.kulhandjian\Desktop\Projects\RadarProject\Code\Matlab\radarDetection\DrondeData\data\Oscar\Results2024';
    outputPathDir = pwd;
end

switch runType
    case "GENDATASET"
        start = clock;
        fprintf( 'Starting Generating Radar Dataset \n');
        % It generates radar dataset according to radar parameters, einvironment.
        generate_radar_dataset(inputPathDir, outputPathDir, mode, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Generating Radar Dataset %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    case "DETECTION"
        start = clock;
        fprintf( 'Starting processRadarDetection \n');
        % It takes the SEC data and processes radar detection
        processRadarDetection( [], inputPathDir, outputPathDir, mode, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished processRadarDetection %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    case "ROC"
        start = clock;
        fprintf( 'Starting ROC Analysis \n');
        % It analyzes Radar Detection via ROC plots.
        roc_analysis( inputPathDir, outputPathDir, mode, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished ROC Analysis %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    case "ESTIMATION"
        start = clock;
        fprintf( 'Starting processRadarEstimation \n');
        % It takes the IQ data and processes radar parameter and CSI estimations
        processRadarEstimation( [], inputPathDir, outputPathDir, mode, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished processRadarEstimation %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    case "ESTIMATION_ANALYSIS"
    otherwise
        disp("Option not available");
end


end

