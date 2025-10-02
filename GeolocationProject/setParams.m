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
	spectrogramFolder   = [resultsFolder '\spectogram' ];
    nameGeolocateRes    = [resultsFolder '\geolocateResult.mat'];
    resultsFile         = [resultsFolder '\results.txt'];
else
    resultsFolder       = [outputPathDir '/results'];
	spectrogramFolder   = [resultsFolder '/spectogram' ];
    nameGeolocateRes    = [resultsFolder '/geolocateResult.mat'];
    resultsFile         = [resultsFolder '/results.txt'];
end

if ~exist(resultsFolder, 'dir')
    % Folder does not exist so create it.
    mkdir(resultsFolder);
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
rP.spectrogramFolder = spectrogramFolder;
rP.nameGeolocateRes  = nameGeolocateRes;
rP.resultsFile       = resultsFile;


rP.isSingleFile      = isSingleFile;
rP.isSpectrogramF    = isSpectrogramF;
rP.isDebugImages     = isDebugImages;
rP.checkVisually     = checkVisually;
rP.saveSigToMat      = saveSigToMat;
rP.isIQtoSECMatlab   = isIQtoSECMatlab;


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
data.sig             = [];
%data.noise           = zeros(Ns,1);

% System Parameters
rP.fs                = fs;
rP.fc                = fc;
rP.bw                = bw;
rP.dataType          = dataType;
rP.dataScalingIQ     = dataScalingIQ;
rP.dataScalingSEC    = dataScalingSEC;
rP.Ns                = Ns;
rP.timeVec           = [0:1:N-1]/fs;
rP.channels          = [];

 
rP.fcLoc            = fcLoc;     % Local Center frequency (MHz) e.g., -15 MHz
rP.pwr              = pwr ;      % Power  Peak power (dB) e.g., -40 dB
rP.beta             = ebeta ;     % Excess Bandwidth e.g.,   0.05
rP.NFFT             = NFFT ;     % Size of FFT  e.g., -32768
rP.out_cmplx        = out_cmplx; % Output IQ e.g., CMPLX' or 'REAL'


rP.pPos              = pPos;
rP.pVel              = pVel;
rP.pUEPos            = pUEPos;
rP.pUEVel            = pUEVel;
rP.bPos              = bPos;
rP.bVel              = bVel;

% Tx Antenna Type
rP.txAntType         = txAntType;

% Tx Antenna Array Params
rP.txArrayDim        = txArrayDim;     

% Rx Antenna Type
rP.rxAntType         = rxAntType;

% Rx Antenna Array Params
rP.rxArrayDim        = rxArrayDim;

%% Channel model
rP.chanModel         = chanModel;

%% add extra AWGN noise based on SNR below
rP.isAWGNChannel     = isAWGNChannel;

%% Generate Radar Dataset parameters
% Parameters
rP.Ts                = Ts; % length of the radar stream in sec
rP.SNR               = SNR;
rP.radioCh           = radioCh;  % Indeces of the rx radios
rP.WAVEFORMS         = WAVEFORMS;












