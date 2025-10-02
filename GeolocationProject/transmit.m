function [rxSig_] = transmit (sig, chanModel, fs, fc, rP, isChannel)

if ~exist('isChannel','var')
    isChannel = true;
end

pPos= rP.pPos;  pVel =rP.pVel;  pUEPos =rP.pUEPos;  pUEVel =rP.pUEVel;
%% Create BS Platform Objects
bPos = rP.bPos;
bVel = rP.bVel;

if length(pPos )<3
    pPos = [pPos ; 0];
    pVel =[pVel ; 0];
    bPos = [bPos; zeros(1,size(bPos,2))];
    bVel = [bVel; zeros(1,size(bVel,2))];
end

%% Create the Tx System platform
[transmitter, radiator, collector, receiver, channel] = setupTxSystem(chanModel, fs, fc, rP);

is5G  = false; % if 5G transmission enabled

if is5G
    %% Create 5G waveform
    [txDL] = setup5GWaveform( fc);
    lSig = length(sig);
    lTxDL = length(txDL);
    if lSig > lTxDL
        txDL = repmat(txDL, ceil(lSig/lTxDL), 1);
    end
    txDL = txDL(1:lSig);
    %txDL = txDL/norm(txDL)*norm(sig);
    txDL = txDL*4;
end

numRx = size(bPos,2 );
[txsig, txstatus] = transmitter(sig);        % Transmit pulse

for idxRx = 1:numRx
%idxRx =1;
    % Calculate the target angles as seen by the BS
    if (isequal(chanModel,'freespace') ||  isequal(chanModel,'two-ray'))
        [bsrngTx, bsangTx] = rangeangle( bPos(:,idxRx) , pPos, chanModel); % Angle between plane antenna and base station antenna
        [bsrngRx, bsangRx] = rangeangle( pPos, bPos(:,idxRx), chanModel); % Angle between plane antenna and base station antenna

        if is5G
            % 5G transmit and receive directions
            [UErngTx,UEangTx]       = rangeangle(bPos(:,idxRx), pUEPos, chanModel);  % second input is the reference
            [UErngRx,UEangRx]       = rangeangle(pUEPos, bPos(:,idxRx), chanModel);
        end
    else
        [bsrngTx, bsangTx] = rangeangle( bPos(:,idxRx), pPos);
        [bsrngRx, bsangRx] = rangeangle( pPos, bPos(:,idxRx));

        if is5G
            % 5G transmit and receive directions
            [UErngTx,UEangTx]       = rangeangle(bPos(:,idxRx), pUEPos);  % second input is the reference
            [UErngRx,UEangRx]       = rangeangle(pUEPos, bPos(:,idxRx));
        end
    end
    % (2) Simulate propagation of pulse in direction of BS
    
    txsigrad             = radiator(txsig, bsangTx);  % Radiate signal

    if is5G
        % UE
        [txsigUE,txstatusUE]  = transmitter(txDL);
        txsigUE               = radiator(txsigUE,UEangTx);
    end

    % (3) Propagate signal to BS
    if isChannel
        % Apply channel model
        if (isequal(chanModel,'freespace') ||  isequal(chanModel,'two-ray'))
            txsigch         = channel( txsigrad, pPos, bPos(:,idxRx), pVel, bVel(:,idxRx));
            if is5G
                txsigUE        = channel(txsigUE,pUEPos,bPos,pUEVel,bVel(:,idxRx));
            end
        else
            [txsigch,~]     = channel( txsigrad);
            if is5G
                [txsigUE,~]   = channel(txsigUE);
            end
        end
    end
    if ~is5G
        rxcol             = collector( txsigch, bsangRx);% Collect
    else
        rxcol          = collector([txsigch txsigUE],[bsangRx UEangRx]);  % this is to receive Signal with 5G
    end
    rxSig             = receiver( rxcol);        % Receiver pre-amp
    %reset(receiver);
    % adding AWGN noise
    %rxSig = awgn(rxSig,awgndB,'measured');

    % if size(rxSig,2) > length(radioCh)
    %     rxSig_ = rxSig(:,radioCh);
    % else
    %     rxSig_ = rxSig;
    % end
    rxSig_(:,:,idxRx) = rxSig;
end

end

function [transmitter, radiator, collector, receiver, channel] = setupTxSystem(chanModel, fs, fc, rP)

c               = physconst('LightSpeed');
lambda           = c/fc;

%% Set UP Transmitter
% Tx Antenna Type
if isfield(rP,'txAntType')
    txAntType = rP.txAntType;
else
    txAntType = 'isotropic'; % Options: {'isotropic', 'URA', 'ULA',  }
end

PLOT_TX_ANT   = false;

% Tx Antenna Array Params
if isfield(rP,'txArrayDim')
    txArrayDim = rP.txArrayDim;
else
    txArrayDim = [1, 1];     % Shape
end

% Setup Tx Antenna
txAntenna     = setupAntenna(txAntType, txArrayDim, fc, lambda, PLOT_TX_ANT);

%% Set Up Transmitter
% Set Up Radiator System Object Properties
peakpower    = 10;   % 10 Watts
txgain       = 10.0;    % Transmitter gain 36 dB (20)

transmitter   = phased.Transmitter( 'PeakPower',peakpower, 'Gain',txgain, 'InUseOutputPort',true);
%% Radiator
radiator     = phased.Radiator('Sensor', txAntenna,'OperatingFrequency', fc,'CombineRadiatedSignals', true);

% Rx Antenna Type
if isfield(rP,'rxAntType')
    rxAntType = rP.rxAntType;
else
    rxAntType = 'isotropic'; % Options: {'isotropic', 'URA', 'ULA', '5GNR', '5GNR_Skylark' (only for 2022a+)}
end

% Rx Antenna Array Params
if isfield(rP,'rxArrayDim')
    rxArrayDim = rP.rxArrayDim;
else
    rxArrayDim = [2, 2];      % Shape
end
PLOT_RX_ANT   = false;

% Setup Rx Antenna
rxAntenna     = setupAntenna(rxAntType, rxArrayDim, fc, lambda, PLOT_RX_ANT);

%% Set Up Receiver
rxgain       = 30.0;     % Receiver gain 42 dB (20)
noisefig     = 5;      % Noise figure of receiver

%% Collector
collector    = phased.Collector( 'Sensor', rxAntenna, 'OperatingFrequency', fc, 'Wavefront','Plane' ); %, 'Polarization', 'Combined');

%% ReceiverPremp
%receiver    = phased.ReceiverPreamp( 'SampleRate', fs, 'Gain', rxgain, 'NoiseFigure', noisefig);
receiver    = phased.Receiver('AddInputNoise',true,'Gain',rxgain,'NoiseFigure',noisefig,'SampleRate',fs);

%% Channel model
grndRFCoeff   = 0.9;
[channel,  nRays] = setupChannel( chanModel, grndRFCoeff, fs, fc);

end

% Setup the Tx antenna and gNB receiver antenna
function antenna = setupAntenna( antType, arrayDim, fc, lambda, view)

% Initialization
frange = [2/3*fc 4/3*fc];

% Create TX/RX Antennas
switch antType
    case 'isotropic'
        antenna = phased.IsotropicAntennaElement;
    case 'ULA'
        % configure isotropic element as the building block
        ant = phased.IsotropicAntennaElement('BackBaffled',false);
        antenna = phased.ULA('Element',ant,'NumElements',max(arrayDim),'ElementSpacing',lambda/2,'ArrayAxis','x');
    case 'URA'
        % Uniform rectangular array (URA)
        antenna = phased.URA('Size', arrayDim, 'ElementSpacing',[lambda/2 lambda/2], 'ArrayNormal','x');
        % Array Pattern
        if view,
            figure;
            pattern(antenna,fc,[-180:180],[-90:90],...
                'CoordinateSystem','polar',...
                'Type','directivity');
            title('Array Pattern');
        end
    case 'URA_A'
        % This example shows how to model a pyramidal conformal antenna array made
        % of 4 panels.
        %% Replicated Array
        sSubA = phased.URA([4,4],'ElementSpacing',[.05 .05]); % Subarray
        set(sSubA.Element,'FrequencyRange',frange,'BackBaffled', true);
        faceSlope = 45;
        antenna = phased.ReplicatedSubarray('Subarray',sSubA,'GridSize',[2 2],...
            'Layout','custom','SubarrayPosition',[.2 0 -.2 0;0 .2 0 -.2; 0 0 0 0],...
            'SubarrayNormal',[0 90 180 -90;90-faceSlope 90-faceSlope 90-faceSlope 90-faceSlope]);
        if view,
            atitle = '45-deg Pyramidal Array - Isotropic Antennas';
            hax = plotArray(antenna,atitle,[]);
        end

        %% Cosine antennas
        antenna.Subarray.Element = phased.CosineAntennaElement(...
            'FrequencyRange',frange,...
            'CosinePower',[8 8]);
        if view,
            atitle = '45-deg Pyramidal Array - Cosine Antennas';
            hax = plotArray(antenna,atitle,hax);
        end

        %% Change face slope
        faceSlope = 75;
        sAnt.SubarrayNormal = [0 90 180 -90;90-faceSlope 90-faceSlope 90-faceSlope 90-faceSlope];
        if view,
            atitle = '75-deg Pyramidal Array - Cosine Antennas';
            hax = plotArray(antenna,atitle,hax);
        end
    case '5GNR'
        antenna = phased.NRAntennaElement('FrequencyRange', [3.1e9 3.45e9], 'Beamwidth', [120 120], 'MaximumGain', 9);
    case '5GNR_Skylark'
        antenna = phased.NRRectangularPanelArray('Size',[1, 7, 1, 6], ...
            'Spacing',[0.5*lambda,0.5*lambda,6*lambda,6*lambda]);
        if view,
            figure; pattern(array,fc,'Orientation',[80;30;60],'ShowArray',true);
        end
    otherwise
        error('Antenna Type Not Supported');
end

end

% Setup Channel
function [channel,  nRays] = setupChannel( chanModel, grndRFCoeff, fs, fc)

%% Free-Space Channel in Two-Way Propgation Model
switch chanModel
    case 'freespace'
        channel = phased.FreeSpace( ...
            'SampleRate',fs, ...
            'OperatingFrequency',fc, ...
            'TwoWayPropagation',false); % One way transmission only
        nRays = 1;
    case 'two-ray'
        channel = twoRayChannel( ...
            'SampleRate',fs, ...
            'OperatingFrequency',fc, ...
            'GroundReflectionCoefficient',grndRFCoeff, ...
            'EnablePolarization',false, ...
            'CombinedRaysOutput',false);
        nRays = 2;
    case 'rician'
        pathDelays = [1.167, 1.258, 1.56, 1.755, 2.987]*1e-6;  %[0 200 800 1200 2300 3700]*1e-9; % in seconds
        avgPathGains = [-3, -9, -10, -23.9, -50]; %[0 -0.9 -4.9 -8 -7.8 -23.9];   % Average path gains (dB)
        KFactor = 3;                   % Linear ratio of specular to diffuse power
        specDopplerShift = 1.6e+03;    % Doppler shift of specular component (Hz)
        fD = 2e+03;                      % Max Doppler shift in Hz
        channel = comm.RicianChannel( ...
            'SampleRate',fs, ...
            'PathDelays',pathDelays, ...
            'AveragePathGains',avgPathGains , ...
            'KFactor', KFactor, ...
            'DirectPathDopplerShift', specDopplerShift, ...
            'MaximumDopplerShift',fD, ...
            'PathGainsOutputPort',true);
        nRays = length(pathDelays);

    otherwise
        error('Channel Model Not Supported');
end

end