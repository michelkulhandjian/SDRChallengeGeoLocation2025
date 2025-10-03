function [data, rP] = setParams( NameFile, inputPathDir, outputPathDir, rP)

% loads the parameters (this part could be a .txt file, it was easier to
% create .m file
params

% Setup the folders
isSingleFile = false;
if exist('NameFile','var')
    % Second parameter exist
    if ~isempty(NameFile)
        isSingleFile = true;
    end
end

if exist('inputPathDir','var')
    % third parameter does exist
    if isempty(inputPathDir)
        inputPathDir = [pwd];
    elseif ~exist(inputPathDir,'dir')
        msg = ['The input directory: ' inputPathDir ',   not found, please verify the directory'];
        error(msg)
    end
else
    inputPathDir = [pwd];
end

if exist('outputPathDir','var')
    % fourth parameter does exist
    if isempty(outputPathDir)
        outputPathDir = [pwd];
    elseif ~exist(outputPathDir,'dir')
        mkdir(outputPathDir);
    end
else
    outputPathDir = [pwd];
end


% Output folders
if ispc
    resultsFolder       = [outputPathDir '\results'];
    secFolder           = [outputPathDir '\SEC'];
    imageFolder         = [resultsFolder '\dataCapture'];
    spectrogramFolder   = [resultsFolder '\spectogram' ];
    nameProcSECData     = [resultsFolder '\proc_sec_data.mat'];
	nameProcEstData     = [resultsFolder '\proc_est_data.mat'];
    resultsFile         = [resultsFolder '\results.txt'];
    nameDetRes          = [resultsFolder '\det_res.mat'];
    nameEstRes          = [resultsFolder '\est_res.mat'];
    nameRCSIJson        = [resultsFolder '\csi_wf_as_pilot.json'];
else
    resultsFolder       = [outputPathDir '/results'];
    secFolder           = [outputPathDir '/SEC'];
    imageFolder         = [resultsFolder '/dataCapture'];
    spectrogramFolder   = [resultsFolder '/spectogram' ];
    nameProcSECData     = [resultsFolder '/proc_sec_data.mat'];
	nameProcEstData     = [resultsFolder '/proc_est_data.mat'];
    resultsFile         = [resultsFolder '/results.txt'];
    nameDetRes          = [resultsFolder '/det_res.mat'];
    nameEstRes          = [resultsFolder '/est_res.mat'];
    nameRCSIJson        = [resultsFolder '/csi_wf_as_pilot.json'];
end

if ~exist(resultsFolder, 'dir')
    % Folder does not exist so create it.
    mkdir(resultsFolder);
end

if ~exist(secFolder, 'dir')
    % Folder does not exist so create it.
    mkdir(secFolder);
end

if ~exist(imageFolder, 'dir')
    % Folder does not exist so create it.
    mkdir(imageFolder);
end

if ~exist(spectrogramFolder, 'dir') && isSpectrogramF
    % Folder does not exist so create it.
    mkdir(spectrogramFolder);
end

% Open a results text file
fileID          = fopen(resultsFile ,'a');
rP.fileID       = fileID;


% Directories
rP.inputPathDir      = inputPathDir;
rP.outputPathDir     = outputPathDir;
rP.resultsFolder     = resultsFolder;
rP.secFolder         = secFolder;
rP.imageFolder       = imageFolder;
rP.spectrogramFolder = spectrogramFolder;
rP.nameProcSECData   = nameProcSECData;
rP.nameProcEstData   = nameProcEstData;
rP.resultsFile       = resultsFile;
rP.nameDetRes        = nameDetRes;
rP.nameEstRes        = nameEstRes;
rP.nameRCSIJson      = nameRCSIJson;

rP.isSingleFile      = isSingleFile;
rP.isSpectrogramF    = isSpectrogramF;
rP.isDebugImages     = isDebugImages;
rP.checkVisually     = checkVisually;
rP.saveSigToMat      = saveSigToMat;
rP.isIQtoSECMatlab   = isIQtoSECMatlab;
rP.maxScheme         = maxScheme;  % 0-original, 1-antenna level, 2-antenna but with location

% Those parameter are set from upper layer
if isfield(rP,'nAntsperGrp')
    nAntsperGrp = rP.nAntsperGrp;
end

if isfield(rP,'WINDOW')
    WINDOW = rP.WINDOW;
end

% parameter initilizaitons (given)
data.totalAntennas   = L*K;
data.nAntsperGrp     = nAntsperGrp;
data.nFreqsperGrp    = nFreqsperGrp;
data.WINDOW          = WINDOW;
data.rWIN            = WINDOW;

% Setting antenna grouping
[data] = setAntGrouping( data);

% store data
data.D               = D;
data.N               = N;
data.fs_SEC          = fs_SEC;
data.num_sparse_elem = num_sparse_elem;

% Initialize Variables
data.numS            = 0;  % Number of Signal Detected
data.numNoise        = 0;  % Number of Noise Detected
data.numSignals      = 0;  % Number of total signals
data.sigInd          = [];  % Signal Index
data.NoiseInd        = [];  % Noise Index
data.sigChannel      = [];  % Signal Index
data.NoiseChannel    = [];  % Noise Index
data.sig             = zeros(Ns,1);
data.noise           = zeros(Ns,1);

% Radar Parameters
rP.fs                = fs;
rP.fc                = fc;
rP.searchGOnly       = true;
rP.antInteg          = true;
rP.dataType          = dataType;
rP.dataScalingIQ     = dataScalingIQ;
rP.dataScalingSEC    = dataScalingSEC;
rP.threshold         = threshold; %0.0025; %2e-5;  0.0059; 0.0055 best range (0.0050 - 0.0065)
rP.DetThreshold      = DetThreshold;%17;  % threshold for radar detection
rP.margin            = margin1;
rP.minAmp            = minAmp; %0.006; %4e-5;
rP.Ns                = Ns;
rP.timeVec           = [0:1:N-1]/fs;
rP.DetInd            = [];
rP.channels          = [];

rP.pPos              = pPos;
rP.pVel              = pVel;
%% Channel model
rP.chanModel         = chanModel;

%% Generate Radar Dataset parameters
% Parameters
rP.Ts                = Ts; % length of the radar stream in sec
rP.SNR               = SNR;
rP.radioCh           = radioCh;  % Indeces of the rx radios
rP.WAVEFORMS         = WAVEFORMS;

%% Radar Detection parameters
rP.DETSNR           = DETSNR;
rP.DETWaveform      = DETWaveform;

%% ROC Analysis parameters
rP.ROCsnr            = ROCsnr;
rP.ROCnAntsperGrp    = ROCnAntsperGrp;
rP.ROCWindow         = ROCWindow;
rP.ROCNspelem        = ROCNspelem;
rP.ROCWaveform       = ROCWaveform;
rP.ROCFig            = ROCFig;










