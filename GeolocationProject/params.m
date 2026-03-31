
%% Common Radar paramaters

% Parameters required for Radar Processing
isSpectrogramF   = false;
isDebugImages    = false;
checkVisually    = false;
saveSigToMat     = false;
isIQtoSECMatlab  = true;

% System parameter initilizaitons
L               = 14; % Number of Antennas per RU mMIMO
K               = 6;  % Number of RUs mMIMO
nAntsperGrp     = 84; % Number of Antennas per Group
nFreqsperGrp    = 84;  % Constrain is 4 Freq per Ant
fs_SEC          = 1;  % sampling rate for SEC
D               = 16;      % Delay Tap
fs              = 100e6; %61.44e6; % Rx sampling rate;          % e.g., 40 MHz
fc              = 3.65e9; %3.6e9;   % center frequency          % e.g., 2.5 GHz 
bw              = 5e6;      % Bandwidth (Hz) of the signal    % e.g., 10 MHz
num_sparse_elem = 1;       % sparseness for SEC
WINDOW          = [30];    % Reporting WIndow size
%numW            = 10;     % number of reporting windows (number of MFs)
N               = [3072];  % Number samples for longest waveform
Ns              = 30720; %2457600; %2457600/4; %6000;    % Max size of samples to be processed by Detectors

%% Modulation Transmitter Parameters 
%params
fcLoc = 0;            % Local Center frequency (MHz) e.g., -15 MHz
pwr = -40;            % Power  Peak power e.g., -40 dB
ebeta = 0.05;         % Excess Bandwidth e.g.,   0.05
NFFT =  32768;        % Size of FFT  e.g., -32768
out_cmplx = 'CMPLX';  % Output IQ e.g., CMPLX' or 'REAL'


dataType        = 'int16';  % Set the data read type (e.g., 'int16', 'float32', etc.)
dataScalingIQ   = 2^(-15);  % This is scaling for IQ data
dataScalingSEC  = 2^(-11);  % This is for SEC output to adhere 5.11
		   
bPos = [0 1.27 2.6; ...
        0 0  1.2]*50; 
	

% bPos = [0 50 75  100 ; ...
%         25 0 30 25 ; ...
%         -52 -60 0 11] ;

bVel = zeros(size(bPos));
%imax = 1700;
%step=100;
%num = 8;
%rng('shuffle');
%bPos = randi([-imax,imax], 3,num);
targetLoc       = 1; % 1- Target Location , 2- Drone, 3- AWACS

if targetLoc   == 1
pPos      = [0.8382; -1.084]*50; %[80; 40; 110];  %Dist 10,409 meters    % 9150 m      % position (m);   At 30,000 ft
pVel      = [0; 0;];  % 
%pPos      = [-20; 20; 10];
%pVel      = [0; 0; 0];  % 
pUEPos    = [15; 20; 10];
pUEVel     = [0; 0; 0];
elseif targetLoc   == 2 % select Drone
   % AWACS: assume aircraft is moving with constant velocity
    pPos      = [150; 100; 300];  %Dist 10,409 meters    % 9150 m      % position (m);   At 30,000 ft
    pVel      = [0; 0;-160];  % -160.934    % (m/s); Travels at 360 mph
    pUEPos    = [15; 20; 10];
    pUEVel     = [0; 0; 0];
elseif targetLoc   == 3 %select AWACS system
    % AWACS: assume aircraft is moving with constant velocity
    pPos      = [1500; 1000; 10e3];  %Dist 10,409 meters    % 9150 m      % position (m);   At 30,000 ft
    pVel      = [113; 113;-20];  % -160.934    % (m/s); Travels at 360 mph
    pUEPos    = [25; 25; 15];
    pUEVel     = [0; 0; 0];
end

% Tx Antenna Type
txAntType     = 'isotropic'; % Options: {'isotropic', 'URA', 'ULA',  }

% Tx Antenna Array Params
txArrayDim    = [1, 1];     % Shape


% Rx Antenna Type
rxAntType     = 'isotropic'; % Options: {'isotropic', 'URA', 'ULA', '5GNR', '5GNR_Skylark' (only for 2022a+)}

% Rx Antenna Array Params
rxArrayDim    = [2, 2];      % Shape

%% Channel model
chanModel     = 'freespace';  % Options: {'freespace', 'two-ray', 'rician', }

%% add extra AWGN noise based on SNR below
isAWGNChannel = false;

%% Generate Radar Dataset Parameters
Ts              = 20e-3; %1e-4 ;%20e-3; % length of the radar stream in sec %  e.g.,  0.1024
SNR             = [30] ; %[-30, -25, -20, -15, -10, 0, 10, 15, 20, 25, 30];
radioCh         = [1:84];  % Indeces of the rx radios
WAVEFORMS       = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B'};
WAVEFORMS      = {'QAM16'}; % {'PSW', 'BPSK', 'QPSK', '8PSK', 'D8PSK', '16PSK', 'QAM16', 'QAM32', 'QAM64', 'QAM256', 'NXDN48', 'NXDN96', 'TONE', 'PAM2', 'PAM4', 'PAM8', 'PAM16'}

%% USRP Hardware Impairment Parameters
%  TX: NI USRP B210 (Analog Devices AD9361 transceiver)
%  RX: NI USRP X310 (UBX-160 daughterboard + ADS62P48 ADC)
usrp.enabled = true;  % Set to false to bypass all USRP impairments

% --- B210 TX Impairments (AD9361 datasheet) ---
usrp.tx.dacBits             = 12;      % AD9361 DAC resolution
usrp.tx.iqGainImbalance_dB  = 0.05;    % After TX quad cal: ~-50 dB image
usrp.tx.iqPhaseImbalance_deg = 0.3;    % After TX quad cal: ~-50 dB image
usrp.tx.dcOffset_dBc        = -50;     % LO leakage after calibration at 0 dB atten
usrp.tx.phaseNoiseRMS_deg   = 0.5;     % Integrated phase noise (~0.37 at 2.4G, ~0.59 at 5.5G)
usrp.tx.enablePA            = true;    % Enable PA nonlinearity model
usrp.tx.oip3_dBm            = 18;      % OIP3 (~19 dBm at 2.4G, ~17 at 5.5G, interpolated for 3.65G)
usrp.tx.outputPower_dBm     = 7;       % Typical TX output power at 3.65 GHz
usrp.tx.maxBW_MHz           = 56;      % Max analog bandwidth (AD9361)

% --- X310 RX Impairments (UBX-160 + ADS62P48) ---
usrp.rx.adcBits             = 14;      % ADS62P48 ADC resolution
usrp.rx.iqGainImbalance_dB  = 0.03;    % UBX-160 RX imbalance (after cal)
usrp.rx.iqPhaseImbalance_deg = 0.2;    % UBX-160 RX phase error
usrp.rx.dcOffset_dBc        = -55;     % RX DC offset after calibration
usrp.rx.phaseNoiseRMS_deg   = 0.4;     % X310 LO phase noise (TCXO reference)
usrp.rx.noiseFigure_dB      = 3;       % UBX-160 noise figure at max gain
usrp.rx.gain_dB             = 31.5;    % UBX-160 max RX gain

% --- Asynchronous Receivers ---
% When true, each X310 receiver has independent clock/LO (no shared 10 MHz/PPS).
% When false, receivers are perfectly synchronized (shared reference).
usrp.rx.asyncEnabled        = true;
usrp.rx.clockOffset_ppm     = 2.5;    % Max sampling clock offset per RX (X310 TCXO: +/-2.5 ppm)
usrp.rx.freqOffset_Hz       = 100;    % Max carrier frequency offset per RX (Hz)
usrp.rx.maxTimeOffset_samples = 50;   % Max random start-time offset (samples) between receivers
usrp.rx.randomPhaseOffset   = true;   % Each RX LO starts at random phase [0, 2*pi)




