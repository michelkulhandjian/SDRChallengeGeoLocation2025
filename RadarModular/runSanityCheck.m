function runSanityCheck (inputPathDir, outputPathDir, sourcePathDir, runType)

%close all; clear all;

% Enable debug logs
rP.debugLog       = true;

if ~exist('runType','var') || (exist('runType', 'var') && isempty(runType))
    runType   = "DETECTION";    % Options: { "SEC", "DETECTION",  "ESTIMATION" }
end

if rP.debugLog
    fprintf( 'Running SanityCheck as runType: %s  \n', runType);
end

if ~exist('inputPathDir','var') || (exist('inputPathDir', 'var') && isempty(inputPathDir)) % Directory in Hobos
    inputPathDir = '/data/data_share/RadarGoldenDatasets/BIN1-A';
end

if ~exist ('outputPathDir', 'var') || (exist('outputPathDir', 'var') && isempty(outputPathDir))
    outputPathDir = '/home/michel.kulhandjian/C++/BIN1ASEC';
end
if ~exist(outputPathDir,'dir')
    mkdir(outputPathDir);
end
if ~exist('sourcePathDir','var') || (exist('sourcePathDir', 'var') && isempty(sourcePathDir)) % Directory in Hobos
    sourcePathDir = '/home/michel.kulhandjian/C++/proj-nsc-rice';
end

curOutputSECDir = outputPathDir;
%curOutputSECSparseDir = [outputPathDir, '/SECSparse' ];
%NameSECMatFile = [outputPathDir,  '/radarSEC.mat'];
rP.source = [sourcePathDir, '/code/radar-detector'];
rP.build = [rP.source, '/build'];
rP.exper_file = [sourcePathDir, '/code/radar-detector/experiment.json'];

switch runType
    case "SEC"
        start = clock;
        fprintf( 'Starting SEC Sanity Check \n');

        % SEC dumper C++
        rP.sec_dumper = [rP.build '/sec_dumper'];
        cmdEx   = [rP.sec_dumper, ' --iq_sample_dir ', inputPathDir, ' --sec_sample_dir ', curOutputSECDir, ' --json_config ',  rP.exper_file];
        system(cmdEx);
        %./build/sec_dumper --json_config experiment.json --iq_sample_dir /data/data_share/RadarGoldenDatasets/BIN1-A --sec_sample_dir /data/data_share/RadarGoldenDatasets/BIN1-A</pre>

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished SEC Sanity Check %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    case "DETECTION"
        start = clock;
        fprintf( 'Starting Coarse Detector Sanity Check \n');

        % corseDetect C++
        % Full SEC--------------------------------------------------------
        %cmake .. -DDISABLE_SPARSENESS=True -S /home/michel.kulhandjian/C++/proj-nsc-rice/code/radar-detector -B /home/michel.kulhandjian/C++/proj-nsc-rice/code/radar-detector/build
        cmdFullSECEx   = [' cmake .. -DDISABLE_SPARSENESS=True', ' -S ', rP.source, ' -B ', rP.build ];
        system(cmdFullSECEx);

        %cmake .. -DDISABLE_SPARSENESS=False -S /home/michel.kulhandjian/C++/proj-nsc-rice/code/radar-detector -B /home/michel.kulhandjian/C++/proj-nsc-rice/code/radar-detector/build
        cmdSparseSECEx   = [' cmake .. -DDISABLE_SPARSENESS=False', ' -S ', rP.source, ' -B ', rP.build ];

        %make -j -C /home/michel.kulhandjian/C++/proj-nsc-rice/code/radar-detector/build
        cmdMakeJEx   = [' make -j', ' -C ', rP.build ];
        system(cmdMakeJEx);

        rP.coarse_detect_iq = [rP.build '/coarse_detect_iq'];
        cmdCoarseDetectEx   = [rP.coarse_detect_iq, ' --iq_sample_dir ' inputPathDir, ' --debug_output true ', ' --json_config ',  rP.exper_file];
        %cmdEx   = [rP.coarse_detect_iq, ' --iq_sample_dir ' inputPathDir, ' --sec_sample_dir ', curOutputSECDir, ' --json_config ',  rP.exper_file];
        system(cmdCoarseDetectEx);
        load 'FullSec.txt';
        load 'Ch1_Sec.txt';
        load 'CoarseDetect_ouput.ascii';

        % Matlab Processing
        % Setting parameters
        rP.Np        = 60;  % this is for fixed W for processing
        rP.numWperPW = 1;

        [data, rP]   = setParams([], inputPathDir, outputPathDir, rP);
        [data, rP]   = loadDataset( inputPathDir, data, rP);
        [Ns, NmF]    = size(data.sig);
        Ns           = 30720;  % this can be disabled if we need to run whole 40ms, this is done to save processing time
        rP.indSig    = 1:Ns;
        [rxSigN]     = data.sig(1:Ns,:);
        rP.maxScheme = 0; % 0-original, 1-antenna level, 2-antenna but with location
        % input Matlab rx, run SEC, Max using reporting window (W) and sparsness (SP) at antenna level (params) - raturn matrices (SEC Full, Sparse
        [data, rP]   = preprocessRadar(data, rxSigN, rP, 1);

        % Compare Full SEC output Matlab vs C++
        errFullSECMax = max(max(abs(data.outputMat - FullSec(:,1:Ns))));
        

        rP.ispMSs=true;
        rP.mode          = 'DETECT';   % LFM Waveform from the Radar Parameters
        rP = setRadarRange( rP);
        if rP.debugLog
            fprintf(' Starting Radar Detector... \n');
        end

        % here we make sure that FULL SEC is used for detection.
        data.STFT_sparse = data.outputMat;
        [ STFT_s, det] =  radarDetector ( data, rP);
        % compare output of MF Matlab vs C++
        NpMSs = length(det.pMSs);
        errOutputFullSECMax = max(abs(det.pMSs - CoarseDetect_ouput(1:NpMSs)'));
        

        data.STFT_sparse = FullSec(:,1:Ns);
        [ STFT_s, det] =  radarDetector ( data, rP);
        % compare output of MF Matlab vs C++
        NpMSs = length(det.pMSs);
        errOutputFullSECCppMax = max(abs(det.pMSs - CoarseDetect_ouput(1:NpMSs)'));
        

        % Sparse SEC-----------------------------------------------------
        system(cmdSparseSECEx);
        system(cmdMakeJEx);
        system(cmdCoarseDetectEx);
        load 'FullSec.txt';
        load 'Ch1_Sec.txt';
        load 'CoarseDetect_ouput.ascii';

        rP.maxScheme = 1; % 0-original, 1-antenna level, 2-antenna but with location
        % input Matlab rx, run SEC, Max using reporting window (W) and sparsness (SP) at antenna level (params) - raturn matrices (SEC Full, Sparse
        [data, rP]   = preprocessRadar(data, rxSigN, rP, 1);

        % Compare Sparse SEC output Matlab vs C++
        errSparseSECMax = max(max(abs(data.STFT_sparse - FullSec(:,1:Ns))));
        

        [ STFT_s, det] =  radarDetector ( data, rP);
        % compare output of MF Matlab vs C++
        NpMSs = length(det.pMSs);
        errOutputSparseSECMax = max(abs(det.pMSs - CoarseDetect_ouput(1:NpMSs)'));
        

        data.STFT_sparse = FullSec(:,1:Ns);
        [ STFT_s, det] =  radarDetector ( data, rP);
        % compare output of MF Matlab vs C++
        NpMSs = length(det.pMSs);
        errOutputSparseSECCppMax = max(abs(det.pMSs - CoarseDetect_ouput(1:NpMSs)'));
        

        % Results
        fprintf( 'Max error of Full SEC between Matlab and C++ is                     %d\n', errFullSECMax);
        fprintf( 'Max error pMSS output with FULL SEC between Matlab and C++ is       %d\n', errOutputFullSECMax);
        fprintf( 'Max error pMSS output with FULL SEC C++ between Matlab and C++ is   %d\n', errOutputFullSECCppMax);

        fprintf( 'Max error of Sparse SEC between Matlab and C++ is                   %d\n', errSparseSECMax);
        fprintf( 'Max error pMSS output with Sparse SEC between Matlab and C++ is     %d\n', errOutputSparseSECMax);
        fprintf( 'Max error pMSS output with Sparse SEC C++ between Matlab and C++ is %d\n', errOutputSparseSECCppMax);

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished Coarse Detector Sanity Check %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    case "ESTIMATION"
        start = clock;
        fprintf( 'Starting Estimator Sanity Check \n');

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        fprintf( 'Finished stimator Sanity Check %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));

    otherwise
        disp("Option not available");
end


