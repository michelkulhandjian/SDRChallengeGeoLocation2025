close all; clear all;

% Enable debug logs
rP.debugLog       = true;   % Outputs debug comments
isIQtoSEC         = true;   % This flag enables converting IQ samples to SEC
isNoiseSingleSrc  = false;  % We can either have single source noise (fixed power) or different according to each SNR
isRadarDetection  = true;   % This enables radar detection
isPlotFigs        = false;  % Print the results

% Select waveform to perform detection analysis (ROCs)
%WAVEFORMS     = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B'};
WAVEFORMS     = {'BIN1-A'};

% SEC dumper C++
rP.sec_dumper = '/data/data_share/bin/sec_dumper';
rP.exper_file = 'experiment.json'; %/data/data_share/bin/experiment.json';
rP.act_ants = ['[ 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, ' ...
    '19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 44, 45, 46, 47, 48, 49, 50, 51, 52, ' ...
    '53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71 ]'];

% Results to be saved
nameResData        = '/space/michel/modular/results/res/res_dataNew.mat';

% Datasets for detection analysis
data(1).inputPathDirIQ  = '/space/SKLK_DATA/MS22_FlightTest_data/Sensitivity'; %multiple datasets
data(2).inputPathDirIQ  = '/data/Michel/20230616_162400_parking-lot-location-1';  % single dataset

% Use this noise dataset if single source noise is enabled,
noise(1).inputPathDirIQ  = '/space/SKLK_DATA/MS22_FlightTest_data/Sensitivity/20230719_183704_full-txgain83/raw_data/gnb/noise';
%noise(1).outputPathDirSEC = '/data/Michel/results/20230616_163514_parking-lot-location-1/BIN3B_noise';


% Parameters we are tunning for
nAntsperGrp = [14, 28];
WINDOW = [30 60];

% save it to res
res(1).nAntsperGrp = nAntsperGrp;
res(1).WINDOW      = WINDOW;
res(1).dataset     = data;

startA = clock;
fprintf( 'Starting Parameter Tuning \n');

for dd = 1: length(data)
    inputPathDirIQ  = data(dd).inputPathDirIQ;
    if ispc
        pathparts = strsplit(inputPathDirIQ,'\');
        outputPathDirSEC = ['\space\michel\modular\results\' pathparts{4}  '\' pathparts{7}];
        outputPathDirSECNoise = ['\space\michel\modular\results\' pathparts{4}  '\' pathparts{7} '_noise'];
    else
        pathparts = strsplit(inputPathDirIQ,'/');
        outputPathDirSEC = ['/space/michel/modular/results/' pathparts{4}  '/' pathparts{7}];
        outputPathDirSECNoise = ['/space/michel/modular/results/' pathparts{4}  '/' pathparts{7} '_noise'];
    end

    inputPathDirSEC = outputPathDirSEC;
    outputPathDirDet = outputPathDirSEC;
    outputPathDirEst = outputPathDirSEC;

    for a = 1: length(nAntsperGrp)
        rP.nAntsperGrp = nAntsperGrp(a);
        for w = 1: length(WINDOW)
            rP.WINDOW = WINDOW(w);
            fprintf( 'nAntsperGrp = %d, WINDOW = %d \n', rP.nAntsperGrp, rP.WINDOW);

            if isIQtoSEC
                start = clock;
                fprintf( 'Starting Signal processIQtoSEC \n');
                % It takes the IQ sample data and converts it to SEC output
                processIQtoSEC( [], inputPathDirIQ, outputPathDirSEC, rP);
                time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
                elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
                fprintf( 'Finished Signal processIQtoSEC %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
            end

            if isRadarDetection
                % %This part is not ported yet
                start = clock;
                fprintf( 'Starting Signal processRadarDetection \n');
                % It takes the SEC data and processes radar detection
                processRadarDetection( [], inputPathDirSEC, outputPathDirDet, rP);
                time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
                elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
                fprintf( 'Finished Signal processRadarDetection %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
            end

            if ispc
                nameDetRes        = [outputPathDirDet '\results\det_res.mat'];
            else
                nameDetRes         = [outputPathDirDet '/results/det_res.mat'];
            end
            load(nameDetRes);
            sigpMSs = det.pMSs;

            % Noise

            inputPathDirIQNoise  = noise(1).inputPathDirIQ;
            inputPathDirSECNoise = outputPathDirSECNoise;
            outputPathDirDetNoise = outputPathDirSECNoise;
            outputPathDirEstNoise = outputPathDirSECNoise;
            if isIQtoSEC
                start = clock;
                fprintf( 'Starting Noise processIQtoSEC \n');
                % It takes the IQ sample data and converts it to SEC output
                processIQtoSEC( [], inputPathDirIQNoise, outputPathDirSECNoise, rP);
                time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
                elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
                fprintf( 'Finished Noise processIQtoSEC %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
            end

            if isRadarDetection
                % %This part is not ported yet
                start = clock;
                fprintf( 'Starting Noise processRadarDetection \n');
                % It takes the SEC data and processes radar detection
                processRadarDetection( [], inputPathDirSECNoise, outputPathDirDetNoise, rP);
                time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
                elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
                fprintf( 'Finished Noise processRadarDetection %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
            end

            if ispc
                nameDetRes        = [outputPathDirDetNoise '\results\det_res.mat'];
            else
                nameDetRes         = [outputPathDirDetNoise '/results/det_res.mat'];
            end
            load(nameDetRes);
            noisepMSs = det.pMSs;
            lSigpMSs = length(sigpMSs);
            noisepMSs1 = noisepMSs(1:lSigpMSs);
            res(dd).sigpMSs{a}{w} = sigpMSs;
            res(dd).noisepMSs{a}{w} = noisepMSs;
            [res(dd).pd{a}{w}, res(dd).pfa{a}{w}] = roc_curve(sigpMSs, noisepMSs1);

            % Save results after each datasets
            save (nameResData,  'res');
        end
    end
end

time_s = (clock - startA)*[0 0 24*60^2 60.^[2 1 0]]';
elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
fprintf( 'Finished tuning parameters %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

% Plot if enabled
if isPlotFigs
    plotDetRes(res);
end


