function radarWaveforms (inputPathDir, outputPathDir , WAVEFORMS, runType)

% Enable debug logs
rP.debugLog       = true;

if ~exist('runType','var')
    runType   = "GENWAVEFORMS";    % Options: {"GENWAVEFORMS", "DOPPLERWAVEFORMS", "CMPWAVEFORMS"}
end

if ~exist('WAVEFORMS','var') || (exist('WAVEFORMS', 'var') && isempty(WAVEFORMS))
    WAVEFORMS = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B' };
end

if rP.debugLog
    fprintf( 'Running as runType: %s  \n', runType);
end

if ~exist('inputPathDir','var') ||  (exist('inputPathDir', 'var') && isempty(inputPathDir))
    inputPathDir  = pwd;
end
if ~exist('outputPathDir','var') || (exist('outputPathDir', 'var') && isempty(outputPathDir))
    outputPathDir = pwd;
end

[~, rP] = setParams([], inputPathDir, outputPathDir, rP);

sigScaling        = 0.9;  % Scaling of the signal

switch runType
    case "GENWAVEFORMS"
        start = clock;
        fprintf( 'Start Generating Radar Waveforms \n');
        for w = 1: length(WAVEFORMS)
            waveF = WAVEFORMS{w};
            if ispc
                filename = [outputPathDir '\', waveF, '.mat'];
            else
                filename = [outputPathDir '/', waveF, '.mat'];
            end
            [sig, psamps, nsamps, Ns] = generateSig (waveF, rP.fs, rP.Ts,  sigScaling);
            save(filename, 'sig', 'nsamps', 'psamps', 'sigStruct');
        end

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Generating Radar Waveforms %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    case "DOPPLERWAVEFORMS"
        start = clock;
        fprintf( 'Start Generating Doppler-Shifted Radar Waveforms \n');

        isMoveAway = false;  % is AWACS moving away (true) or coming towards (false) ?
        % 1) AWACS - about 1.8kHz max doppler shift under the following conditions.
        % 2) Drone - about 70 Hz max doppler shift under the following conditions.
        fd_a    = 1800;
        fd_d    = 70;
        fD      = fd_a - fd_d;    % Doppler Frequency

        if isMoveAway
            fD = -fD;
        end

        for w = 1: length(WAVEFORMS)
            waveF = WAVEFORMS{w};
            if ispc
                inputfilename = [inputPathDir '\', waveF, '.mat'];
                filename = [outputPathDir '\', waveF, '.mat'];
            else
                inputfilename = [inputPathDir '/', waveF, '.mat'];
                filename = [outputPathDir '/', waveF, '.mat'];
            end

            %[sig, psamps, nsamps, Ns] = generateSig (waveF, rP.fs, rP.Ts,  sigScaling);
            load(inputfilename);
            % Doppler Shifting
            N              = length(sig);
            nRange         = [0:N-1]';
            sig = bsxfun(@times,sig,exp(1i*2*pi*k*fD*(Ts*(nRange))));
            nsamps         = psamps;
            save(filename, 'sig', 'nsamps', 'psamps', 'sigStruct');
        end

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Generating Doppler-Shifted Radar Waveforms %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    case "CMPWAVEFORMS"
        start = clock;
        fprintf( 'Starting Comparing Radar Waveforms \n');
        for w = 1: length(WAVEFORMS)
            waveF = WAVEFORMS{w};
            if ispc
                filename = [inputPathDir '\', waveF, '.mat'];
            else
                filename = [inputPathDir '/', waveF, '.mat'];
            end

            [sig_, Np, PRI, Ns] = generateSig (waveF, rP.fs, rP.Ts,  sigScaling);
            load ( filename) ;
            if ( norm(sig - sig_) ~= 0 || PRI ~=nsamps || psamps ~= Np )
                fprintf( 'Waveform: %s  is not LFM waveform\n', filename);
            else
                fprintf( 'Waveform: %s  is LFM waveform\n', filename);
            end

        end

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Comparing Radar Waveforms %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));


    otherwise
        disp("Option not available");
end