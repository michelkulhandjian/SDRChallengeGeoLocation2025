function [rxSig_] = transmitRadar (sig, chanModel, fs, fc, pPos,  pVel, isChannel)

if ~exist('isChannel','var')
    isChannel = true;
end

%% Create the Radar System platform
[bPos, bVel, transmitter, radiator, collector, receiver, channel] = setupRadarSystem(chanModel, fs, fc);

% Calculate the target angles as seen by the BS
if (isequal(chanModel,'freespace') ||  isequal(chanModel,'two-ray'))
    [bsrngTx, bsangTx] = rangeangle( bPos, pPos, chanModel); % Angle between plane antenna and base station antenna
    [bsrngRx, bsangRx] = rangeangle( pPos, bPos, chanModel); % Angle between plane antenna and base station antenna
else
    [bsrngTx, bsangTx] = rangeangle( bPos, pPos);
    [bsrngRx, bsangRx] = rangeangle( pPos, bPos);
end
% (2) Simulate propagation of pulse in direction of BS
[txsig, txstatus] = transmitter(sig);        % Transmit pulse
txsig             = radiator(txsig, bsangTx);  % Radiate signal

% (3) Propagate signal to BS
if isChannel
    % Apply channel model
    if (isequal(chanModel,'freespace') ||  isequal(chanModel,'two-ray'))
        txsig         = channel( txsig, pPos, bPos, pVel, bVel);
    else
        [txsig,~]     = channel( txsig);
    end
end
rxcol             = collector( txsig, bsangRx);% Collect
rxSig             = receiver( rxcol);        % Receiver pre-amp

% if size(rxSig,2) > length(radioCh)
%     rxSig_ = rxSig(:,radioCh);
% else
%     rxSig_ = rxSig;
% end
rxSig_ = rxSig;

end

function [bPos, bVel,  transmitter, radiator, collector, receiver, channel] = setupRadarSystem(chanModel, fs, fc)

c               = physconst('LightSpeed');
lambda           = c/fc;

%% Set UP Transmitter
% Tx Antenna Type
txAntType     = 'isotropic'; % Options: {'isotropic', 'URA', 'ULA',  }
PLOT_TX_ANT   = false;

% Tx Antenna Array Params
txArrayDim    = [1, 1];     % Shape

% Setup Tx Antenna
txAntenna     = setupAntenna(txAntType, txArrayDim, fc, lambda, PLOT_TX_ANT);

%% Set Up Transmitter
% Set Up Radiator System Object Properties
peakpower    = 10;   % 10 Watts
txgain       = 10.0;    % Transmitter gain 36 dB (20)

transmitter   = phased.Transmitter( 'PeakPower',peakpower, 'Gain',txgain, 'InUseOutputPort',true);
%% Radiator
radiator     = phased.Radiator('Sensor', txAntenna,'OperatingFrequency', fc,'CombineRadiatedSignals', true);

%% Create BS Platform Objects
% BS: assume the base station is stationary at 20m above the ground (xy origin)
bPos         = [0; 0; 0];  % position (m)
bVel         = [0; 0; 0];   % (m/s)

% Rx Antenna Type
rxAntType     = '5GNR_Skylark'; % Options: {'isotropic', 'URA', 'ULA', '5GNR', '5GNR_Skylark' (only for 2022a+)}

% Rx Antenna Array Params
rxArrayDim    = [2, 2];      % Shape
PLOT_RX_ANT   = false;

% Setup Rx Antenna
rxAntenna     = setupAntenna(rxAntType, rxArrayDim, fc, lambda, PLOT_RX_ANT);

%% Set Up Receiver
rxgain       = 30.0;     % Receiver gain 42 dB (20)
noisefig     = 5;      % Noise figure of receiver

%% Collector
collector    = phased.Collector( 'Sensor', rxAntenna, 'OperatingFrequency', fc, 'Wavefront','Plane' ); %, 'Polarization', 'Combined');

%% ReceiverPremp
receiver    = phased.ReceiverPreamp( 'SampleRate', fs, 'Gain', rxgain, 'NoiseFigure', noisefig);

%% Channel model
grndRFCoeff   = 0.9;
[channel,  nRays] = setupChannel( chanModel, grndRFCoeff, fs, fc);

end

% Setup the Airborne radar antenna and gNB receiver antenna
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