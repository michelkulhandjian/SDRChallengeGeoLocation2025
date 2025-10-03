function [rparams] = radarParams(mode)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  Generate the Following
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % % 1A - Min PW, min Bc, and mid-range PRF
    % % 1B - Max PW, max Bc, and mid-range PRF
    % % 2A - Min PW, max Bc, and mid-range PRF
    % % 2B - Max PW, min Bc, and mid-range PRF
    % % 3A - Min PW, min Bc, and mid-range PRF
    % % 3B - Max PW, min Bc, and mid-range PRF
    % % OTHER A - Bin 2A (will add Doppler)
    % % OTHER B - Bin 3A (will add Doppler)

    % UPDATED FEB. 2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Common vals
    rparams.pFA                = 0.001;                   % probability of false alarm
    rparams.TX_SCALE           = 1;%db2mag(-13);             % Scale for Tx waveform ([0:1])
    rparams.N_ZPAD_PRE         = 150;                     % Zero-padding prefix for Iris
    rparams.N_ZPAD_POST        = 150;                     % Zero-padding postfix for Iris
    rparams.N_ZPAD_MID         = 50;

    switch mode
    case 'BIN1-A'
        % Min PW, min Bc, and mid-range PRF
        rparams.fs_radar           = 61.44e6;                       % sample frequency
        rparams.pulseWidth         = 1e-6;
        rparams.prf                = 20e3;
        rparams.BW_hps             = 1e6/1e-6;                  % MHz/usec (microsecond)
        rparams.sweepBW            = rparams.pulseWidth * rparams.BW_hps;     % Hz
        rparams.numPulses          = 20;
     
    case 'BIN1-B'
        % Max PW, max Bc, and max PRF
        rparams.fs_radar           = 61.44e6;                       % sample frequency
        rparams.pulseWidth         = 2e-6;
        rparams.prf                = 30e3;
        rparams.BW_hps             = 4e6/1e-6;                  % MHz/usec (microsecond)
        rparams.sweepBW            = rparams.pulseWidth * rparams.BW_hps;     % Hz
        rparams.numPulses          = 20;

    case 'BIN2-A'
        % Min PW, mid Bc, and max PRF
        rparams.fs_radar           = 61.44e6;                      % Ssample frequency     
        rparams.pulseWidth         = 50e-6;
        rparams.prf                = 200;
        rparams.BW_hps             = nan;
        rparams.sweepBW            = 20e6;                      % Hz
        rparams.numPulses          = 20;

    case 'BIN2-B'
        % Max PW, min Bc, and max PRF
        rparams.fs_radar           = 61.44e6;                      % sample frequency
        rparams.pulseWidth         = 100e-6;
        rparams.prf                = 200;
        rparams.BW_hps             = nan;
        rparams.sweepBW            = 30e6;                      % Hz
        rparams.numPulses          = 20;

    case 'BIN3-A'
        % Min PW, min Bc, and mid-range PRF
        rparams.fs_radar           = 61.44e6;                       % sample frequency
        rparams.pulseWidth         = 10e-6;
        rparams.prf                = 100;
        rparams.BW_hps             = nan;
        rparams.sweepBW            = 5e6;                       % Hz
        rparams.numPulses          = 20;

    case 'BIN3-B'
        % Max PW, min Bc, and mid-range PRF  (PRF must be between 30 and 150 Hz)
        rparams.fs_radar           = 61.44e6;                       % sample frequency
        rparams.pulseWidth         = 50e-6;
        rparams.prf                = 150;
        rparams.BW_hps             = nan;                       % MHz
        rparams.sweepBW            = 15e6;                      % Hz
        rparams.numPulses          = 20;

    case 'OTHER-A'
        % Bin 2B ( ADD DOPPLER!! )
        rparams.fs_radar           = 40e6;                      % Ssample frequency     
        rparams.pulseWidth         = 50e-6;
        rparams.prf                = 100;
        rparams.BW_hps             = nan;
        rparams.sweepBW            = 20e6;                      % Hz
        rparams.numPulses          = 20;

    case 'OTHER-B'
        % Bin 3A ( ADD DOPPLER!! )
        rparams.fs_radar           = 40e6;                       % sample frequency
        rparams.pulseWidth         = 10e-6;
        rparams.prf                = 50;
        rparams.BW_hps             = nan;
        rparams.sweepBW            = 5e6;                       % Hz
        rparams.numPulses          = 20;

    otherwise
        error('Invalid Radar Configuration Option');
    end
    %fprintf("CONFIGURED FOR RADAR MODE: %s \n",mode);
    %disp(rparams);
end
