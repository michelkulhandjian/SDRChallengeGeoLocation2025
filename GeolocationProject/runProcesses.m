function runProcesses (inputPathDir, outputPathDir , runType, mode )

%close all; clear all;

% Enable debug logs
rP.debugLog       = true;

if ~exist('runType','var')
    runType   = "GENDATASET";   % Options: {"GENDATASET", "LOCANALYSIS", "LOCALIZE"}
end

if ~exist('mode','var')
    mode = 1;
end

if rP.debugLog
    fprintf( 'Running as runType: %s with mode %d \n', runType, mode);
end

if ~exist('inputPathDir','var') ||  (exist('inputPathDir', 'var') && isempty(inputPathDir))
    inputPathDir  = 'C:\Users\mkulhandjian\OneDrive - Digital Global Systems\Desktop\codes\TDOA\Matlab\DGSTDOA\datasetSNR30\raw_data\gnb\QAM16';
    %inputPathDir  = pwd;
end
if ~exist('outputPathDir','var') || (exist('outputPathDir', 'var') && isempty(outputPathDir))
    outputPathDir  = 'C:\Users\michel.kulhandjian\Desktop\Projects\RadarProject\Code\Matlab\radarDetection\DrondeData\data\Oscar\Results2024';
    outputPathDir = pwd;
end

switch runType
    case "GENDATASET"
        start = clock;
        fprintf( 'Starting Generating Dataset \n');
        % It generates radar dataset according to radar parameters, einvironment.
        mode = 129;
        generate_dataset(inputPathDir, outputPathDir, mode, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Generating Dataset %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    case "LOCANALYSIS"
         start = clock;
         fprintf( 'Starting localization analaysis \n');
         % It takes the SEC data and processes radar detection
         locanalysis( [], inputPathDir, outputPathDir, mode, rP);
         time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
         elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
         fprintf( 'Finished localization analaysis %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    case "LOCALIZE"
         start = clock;
         fprintf( 'Starting localization \n');
         % It analyzes Radar Detection via ROC plots.
         mode = 1;
         geolocate([], inputPathDir, outputPathDir, mode, rP);
         time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
         elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
         fprintf( 'Finished LOCALIZE %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
    otherwise
        disp("Option not available");
end


end

