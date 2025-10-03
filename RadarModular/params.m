
%% Common Radar paramaters

% Parameters required for Radar Processing
isSpectrogramF   = false;
isDebugImages    = false;
checkVisually    = false;
saveSigToMat     = false;
isIQtoSECMatlab  = true;
maxScheme        = 1;  % 0-original, 1-antenna level, 2-antenna but with location

% Radar parameter initilizaitons
L               = 14; % Number of Antennas per RU mMIMO
K               = 6;  % Number of RUs mMIMO
nAntsperGrp     = 84; % Number of Antennas per Group
nFreqsperGrp    = 84;  % Constrain is 4 Freq per Ant
fs_SEC          = 1;  % sampling rate for SEC
D               = 16;      % Delay Tap
fs              = 61.44e6; % Rx sampling rate;
fc              = 3.6e9;   % center frequency
num_sparse_elem = 1;       % sparseness for SEC
WINDOW          = [30];    % Reporting WIndow size
%numW            = 10;     % number of reporting windows (number of MFs)
N               = [3072];  % Number samples for longest waveform
Ns              = 30720; %2457600; %2457600/4; %6000;    % Max size of samples to be processed by Detectors

dataType        = 'int16';  % Set the data read type (e.g., 'int16', 'float32', etc.)
dataScalingIQ   = 2^(-15);  % This is scaling for IQ data
dataScalingSEC  = 2^(-11);  % This is for SEC output to adhere 5.11

% Thresholds
threshold       = 2e-5; %0.0025; %2e-5;  0.0059; 0.0055 best range (0.0050 - 0.0065)
margin1         = 6;
minAmp          = 2e-5; %0.006; %4e-5;
DetThreshold    = 1.4; %3;%17;  % threshold for radar detection

% AWACS: assume aircraft is moving with constant velocity
if true % select Drone
    pPos      = [150; 100; 300];  %Dist 10,409 meters    % 9150 m      % position (m);   At 30,000 ft
    pVel      = [0; 0;-160];  % -160.934    % (m/s); Travels at 360 mph
else %select AWACS system
    % AWACS: assume aircraft is moving with constant velocity
    pPos      = [1500; 1000; 10e3];  %Dist 10,409 meters    % 9150 m      % position (m);   At 30,000 ft
    pVel      = [113; 113;-20];  % -160.934    % (m/s); Travels at 360 mph
end

%% Channel model
chanModel     = 'freespace';  % Options: {'freespace', 'two-ray', 'rician', }

%% Generate Radar Dataset Parameters
Ts              = 20e-3; %1e-4 ;%20e-3; % length of the radar stream in sec
SNR             = 10; %[-30, -25, -20, -15, -10, 0, 10, 15, 20, 25, 30];
radioCh         = [1:84];  % Indeces of the rx radios
WAVEFORMS       = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B'};
WAVEFORMS      = {'BIN1-A'};


%% Radar Detection parameters
DETSNR           = [-15];
DETWaveform      = 'BIN1-A';

%% ROC Analysis parameters
% Scenarios to plot/present
ROCsnr          = -15; %[-20 -15 -10];%[1:10]; %[5, 8, 10, 20]; % Indexes, not actual gains
ROCnAntsperGrp  = 84; %[1 14 64 84];  % Number of antenna integration per Group
ROCWindow       = 30; %[ 15 30 60]; % MF window size %[32 62 150 185 250]; %[ 25 35 50 61]; [ 32, 64, 85, 100, 122]%
ROCNspelem      = 1; %[1 16 30];  % sparseness for SEC (# max values)
ROCWaveform     = 'BIN1-A';
ROCFig          = 1; % Plotting ROC figure type
% 1. fix SNR, W, AntPerGrp, sp           - One ROC plot
% 2. Sweeping SNR (fix W, AntPerGrp, sp) - One ROC plot (several SNRs)
% 3. Sweeping SNR, W, AntPerGrp, fix sp  - 3x3 ROCs (several SNRs each plot)
% 4. Sweeping W (fix SNR, AntPerGrp, sp) - One ROC plot (several Ws)
% 5. Sweeping W, SNR, AntPerGrp, fix sp  - 3x3 ROCs (several Ws each plot) - prev
% 6. Sweeping AntPerGrp (fix SNR, W, sp) - One ROC plot (several AntPerGrps)
% 7. Sweeping AntPerGrp, SNR, W, fix sp  - 3x3 ROCs (several AntPerGrps each plot)
% 8. Sweeping sp (fix SNR, W, AntPerGrp) - One ROC plot (several sps)
% 9. Sweeping sp, SNR, AntPerGrp, fix W - 3x3 ROCs (several sps each plot) - prev

