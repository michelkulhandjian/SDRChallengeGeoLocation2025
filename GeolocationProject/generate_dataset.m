% Generate Dataset,
function [data, rP] = generate_dataset( inputPathDirIQ, outputPathDirIQ, mode, rP)

% For generating signal, channel parameters, snr are taken from params.m
% Generate signal, apply channel + snr (rxsig)
% Run SEC is to run SEC algorithm plus apply reporting window + delta-max sparsness
% Run AntInt is to run antenna integration on SEC according to antenna grouping set on params.m

% 1. Generate rxsig, save sig.mat,sig.dat + run SEC save, SEC.mat, SECSparse.mat, SEC.dat, SECSParse.dat
%  a. generateSig()
%  b. transmit()
%  c. addingAWGN()
%  d. saveMAT() - Sig
%  e. writeToFile() Sig
%  f. Preprocessing() - mode 2
%  g. saveMAT() - SEC Full + Sparse
%  h. writeToFile() - SEC Full + Sparse mode

% 2. load IQ dataset, save sig.mat + run SEC save, SEC.mat, SECSparse.mat, SEC.dat, SECSParse.dat
%  a. loadDataset()
%  d. saveMAT() - Sig
%  d. Preprocessing() - mode 2
%  e. saveMAT() - SEC
%  f. writeToFile() - SEC

% 3. load SEC dataset, run AntInt save outputMat.mat
%  a. loadDataset()
%  b. Preprocessing() - mode 3

% mode composed of bits: |_7_|_6_|_5_|_4_|_3_|_2_|_1_|_0_|
% 0-bit (1)  represents generate signal, rxsig
% 1-bit (2)  represents write signal to file
% 2-bit (4)  represents run SEC + reporting window + delta-max sparsness
% 3-bit (8)  represents save SEC full into files
% 4-bit (16) represents save SEC sparse into files
% 5-bit (32) represents process IQ dataset
% 6-bit (64) represents run antenna integration
% note that SEC full and SEC sparse saved into .mat files are dependent on saveSigToMat
% Masks applied to mode
% [0 0 0 0 0 0 1] (1)  - rxsig save sig.mat
% [0 0 0 0 0 1 1] (3)  - rxsig save sig.mat, sig.dat
% [0 0 0 0 1 0 1] (5)  - rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat
% [0 0 0 0 1 1 1] (7)  - rxsig save sig.mat, sig.dat + run SEC save SEC.mat, SECSparse.mat
% [0 0 0 1 1 0 1] (13) - rxsig save sig.mat + run SEC save SEC.mat, SEC.dat, SECSparse.mat
% [0 0 0 1 1 1 1] (15) - rxsig save sig.mat, sig.dat + run SEC save SEC.mat, SECSparse.mat, SEC.dat
% [0 0 1 0 1 0 1] (21) - rxsig save sig.mat  + run SEC save SEC.mat, SECSparse.mat, SECSparse.dat
% [0 0 1 0 1 1 1] (23) - rxsig save sig.mat, sig.dat + run SEC save SEC.mat, SECSparse.mat, SECSparse.dat
% [0 0 1 1 1 0 1] (29) - rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat, SECSparse.mat, SEC.dat, SECSparse.dat
% [0 0 1 1 1 1 1] (31) - rxsig save sig.mat, sig.dat + run SEC save SEC.mat, SECSparse.mat, SEC.dat, SECSparse.dat
% [0 1 0 0 0 0 0] (32) - load IQ rxsig save sig.mat
% [0 1 0 0 1 0 0] (36) - load IQ rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat
% [0 1 0 1 1 0 0] (44) - load IQ rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat, SEC.dat
% [0 1 1 0 1 0 0] (52) - load IQ rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat, SECSparse.dat
% [0 1 1 1 1 0 0] (60) - load IQ rxsig save sig.mat + run SEC save SEC.mat, SECSparse.mat, , SEC.dat, SECSparse.dat
% [1 0 0 0 0 0 0] (64) - load SEC save SEC.mat + run AntInt save outputMat.mat

% Common parmaters
[data, rP] = setParams([], inputPathDirIQ, outputPathDirIQ, rP);

if rP.debugLog
    fprintf(' Start Dataset Generating Process... \n');
end
% Write to results
fprintf(rP.fileID,' Start Dataset Generating Process... \n');

% Output folders
if ispc
    gnbPart       = ['\raw_data\gnb'];
else
    gnbPart       = ['/raw_data/gnb'];
end

outputPathDir    = rP.outputPathDir;
inputPathDir     = rP.inputPathDir;

if bitand(mode , 128) % Localization
    % Parameters
    fs                = rP.fs;
    fc                = rP.fc;
    bw                = rP.bw;
    Ts                = rP.Ts; % length of the stream in sec
    SNR               = rP.SNR;
    WAVEFORMS         = rP.WAVEFORMS;

    %% Channel model
    chanModel         = rP.chanModel;
    sigScaling        = 0.9;  % Scaling of the signal
    pdBm              = 0; %-113;

    Nsnr   = length(SNR);
    NwaveF = length(WAVEFORMS)

    for indSnr = 1: Nsnr
        curSNR = SNR(indSnr);
        for w = 1 : NwaveF
            %Select the  Waveform
            curWaveF   = WAVEFORMS{w};   % LFM Waveform from the Parameters
            fprintf('SNR: %d dB (%d/%d), Waveform: %s (%d/%d)\n', ...
                curSNR, indSnr, Nsnr, ...
                curWaveF, w, NwaveF );

            % Create the folder for Rx and noise dataset, and .mat files
            if ispc
                curOutputDir = [outputPathDir, '\datasetSNR', int2str(SNR(indSnr)), gnbPart, '\' , curWaveF];
                NameSigMatFile = [curOutputDir,  '\Sig.mat'];
            else
                curOutputDir = [outputPathDir, '/datasetSNR', int2str(SNR(indSnr)), gnbPart, '/' , curWaveF];
                NameSigMatFile = [curOutputDir,  '/Sig.mat'];
            end
            if ~exist(curOutputDir,'dir')
                mkdir(curOutputDir);
            end

            % Generate the  waveform
            %[sig, Np, PRI, Ns]
            [sig, ~, ~, ~] = generateSig (curWaveF, fs, Ts, bw,  sigScaling, rP);
            [rxSig_] = transmit (sig, chanModel, fs, fc, rP, true); rxSig_ = squeeze(rxSig_);
            %pRMS = rms(rxSig_(:))^2;
            %rxSig_ = rxSig_/sqrt(pRMS)*10^((pdBm - 30)/20);
            numRx = size(rxSig_, 2 );
            for idxRx = 1:numRx
                % Add AWGN noise accross all the received antennas
                if rP.isAWGNChannel
                    [rxSigN, Noise] = addingAWGN (rxSig_(:,idxRx), curSNR);
                else
                    [rxSigN]  = rxSig_(:,idxRx);
                end

                data.sig(:,idxRx) = rxSigN;

                % if 1 bit is set we save it to file
                if bitand(mode , 1)
                    if ispc
                        NameFile = [curOutputDir,  '\rx', sprintf('%02d',idxRx), '_samples.dat'];
                    else
                        NameFile = [curOutputDir,  '/rx', sprintf('%02d',idxRx), '_samples.dat'];
                    end

                    if rP.debugLog
                        fprintf(' file name %s \n', NameFile);
                    end

                    % Save received  waveform in .dat file
                    writeToFile (rxSigN.', NameFile, rP.dataType, rP.dataScalingIQ );
                end
            end
            if  rP.saveSigToMat
                save (NameSigMatFile, 'data', 'rP', '-v7.3');
            end
        end
    end

elseif bitand(mode , 1)

    % Parameters
    fs                = rP.fs;
    fc                = rP.fc;
    bw                = rP.bw;
    Ts                = rP.Ts; % length of the stream in sec
    SNR               = rP.SNR;
    radioCh           = rP.radioCh;  % Indeces of the rx radios
    WAVEFORMS         = rP.WAVEFORMS;

    %% Channel model
    chanModel         = rP.chanModel;
    sigScaling        = 0.9;  % Scaling of the signal
    pdBm              = 0; %-113;

    Nsnr   = length(SNR);
    NwaveF = length(WAVEFORMS)

    for indSnr = 1: Nsnr
        curSNR = SNR(indSnr);
        for w = 1 : NwaveF
            %Select the  Waveform
            curWaveF   = WAVEFORMS{w};   % LFM Waveform from the Parameters
            fprintf('SNR: %d dB (%d/%d), Waveform: %s (%d/%d)\n', ...
                curSNR, indSnr, Nsnr, ...
                curWaveF, w, NwaveF );

            % Create the folder for Rx and noise dataset, and .mat files
            if ispc
                curOutputDir = [outputPathDir, '\datasetSNR', int2str(SNR(indSnr)), gnbPart, '\' , curWaveF];
                curOutputDirNoise = [outputPathDir, '\datasetSNR', int2str(SNR(indSnr)), gnbPart, '\noise' ];
                curOutputSECDir = [outputPathDir, '\datasetSNR', int2str(SNR(indSnr)), gnbPart, '\SEC' ];
                curOutputSECSparseDir = [outputPathDir, '\datasetSNR', int2str(SNR(indSnr)), gnbPart, '\SECSparse' ];
                NameSigMatFile = [curOutputDir,  '\Sig.mat'];
                NameSECMatFile = [curOutputDir,  '\SEC.mat'];
            else
                curOutputDir = [outputPathDir, '/datasetSNR', int2str(SNR(indSnr)), gnbPart, '/' , curWaveF];
                curOutputDirNoise = [outputPathDir, '/datasetSNR', int2str(SNR(indSnr)), gnbPart, '/noise' ];
                curOutputSECDir = [outputPathDir, '/datasetSNR', int2str(SNR(indSnr)), gnbPart, '/SEC' ];
                curOutputSECSparseDir = [outputPathDir, '/datasetSNR', int2str(SNR(indSnr)), gnbPart, '/SECSparse' ];
                NameSigMatFile = [curOutputDir,  '/Sig.mat'];
                NameSECMatFile = [curOutputDir,  '/SEC.mat'];
            end
            if ~exist(curOutputDir,'dir')
                mkdir(curOutputDir);
            end
            if ~exist(curOutputDirNoise,'dir')
                mkdir(curOutputDirNoise);
            end

            % Generate the  waveform
            %[sig, Np, PRI, Ns]
            [sig, ~, ~, ~] = generateSig (curWaveF, fs, Ts, bw, sigScaling, rP);
            [rxSig_] = transmit (sig, chanModel, fs, fc, rP, true);
            pRMS = rms(rxSig_(:))^2;
            rxSig_ = rxSig_/sqrt(pRMS)*10^((pdBm - 30)/20);

            % Add AWGN noise accross all the received antennas
            [rxSigN, Noise] = addingAWGN (rxSig_, curSNR);

            if  rP.saveSigToMat
                data.sig = rxSigN;
                save (NameSigMatFile, 'data', 'rP', '-v7.3');
            end
            %=============================================================
            if bitand( mode , 2 ) % Writing to .dat binary files
                for indCh = 1:length (radioCh)
                    curChan = radioCh(indCh);
                    % Name of the file
                    if ispc
                        NameFile = [curOutputDir,  '\ch', sprintf('%02d',curChan), '_samples.dat'];
                        NameFileNoise = [curOutputDirNoise,  '\ch', sprintf('%02d',curChan), '_samples.dat'];
                    else
                        NameFile = [curOutputDir,  '/ch', sprintf('%02d',curChan), '_samples.dat'];
                        NameFileNoise = [curOutputDirNoise,  '/ch', sprintf('%02d',curChan), '_samples.dat'];
                    end
                    if rP.debugLog
                        fprintf(' file name %s \n', NameFile);
                    end

                    % Save received  waveform in .dat file
                    writeToFile (rxSigN(:,indCh).', NameFile, rP.dataType, rP.dataScalingIQ );
                    % Verify
                    %[x2] = parseBinFile (NameFile, dataType, dataScalingIQ);
                    %norm(rxSigN(:,indCh) - x2)
                    if w == 1
                        % Save noise in .dat file
                        writeToFile (Noise(:,indCh).', NameFileNoise, rP.dataType, rP.dataScalingIQ );
                        % Verify
                        %[x2] = parseBinFile (NameFileNoise, dataType, dataScalingIQ);
                        %norm(Noise(:,indCh) - x2)
                    end
                end % End of radio channels
            end % End of mode 2

            if bitand( mode , 4 ) %  run SEC
                % Run SEC + r.Window and sparseness of the received .mat matrix
                [data] = processSEC(rxSigN, curOutputSECDir, curOutputSECSparseDir, NameSECMatFile, data, mode, rP);
            end % End Mode 4

        end % End of waveforms
    end % mode 1
elseif  bitand(mode , 32) % load IQ dataset
    % Create the SEC, SECSparse folders as well as sig.mat, SEC.mat files
    if ispc
        curOutputSECDir = [outputPathDir, '\SEC' ];
        curOutputSECSparseDir = [outputPathDir, '\SECSparse' ];
        NameSigMatFile = [outputPathDir,  '\Sig.mat'];
        NameSECMatFile = [outputPathDir,  '\SEC.mat'];
    else
        curOutputSECDir = [outputPathDir,  '/SEC' ];
        curOutputSECSparseDir = [outputPathDir, '/SECSparse' ];
        NameSigMatFile = [outputPathDir,  '/Sig.mat'];
        NameSECMatFile = [outputPathDir,  '/SEC.mat'];
    end

    myFiles = dir(fullfile(rP.inputPathDir,'**','*.dat')); %gets all mat files in struct
    NmF = length(myFiles);

    if NmF > 0
        data.totalAntennas   = NmF;
        [data] = setAntGrouping( data);
        % We can either process using Matlab or C++
        if ~rP.isIQtoSECMatlab
            % C++ processing
            % Write to json file
            experiment_to_json( data.nAntsperGrp, data.nFreqsperGrp, rP.act_ants, rP.exper_file);
            cmdEx   = [rP.sec_dumper, ' --iq_sample_dir ' inputPathDir, ' --sec_sample_dir ', curOutputSECDir, ' --json_config ',  rP.exper_file];
            system(cmdEx);
            %system('./sec_dumper --iq_sample_dir /data/data_share/20230616_162400_parking-lot-location-1/raw_data/gnb/BIN1-A --sec_sample_dir /data/data_share/20230616_162400_parking-lot-location-1/sec/BIN1-A --json_config experiment.json')
            % Verification
            if (0)
                baseFileName = myFiles(1).name;
                [ext, name] = fileparts(baseFileName);
                if ispc
                    pathparts = strsplit(myFiles(1).folder,'\');
                    NameFile = [myFiles(1).folder '\' myFiles(1).name];
                    NameSECFile = [curOutputSECDir '\' name '.sec'];
                else
                    pathparts = strsplit(myFiles(1).folder,'/');
                    NameFile = [myFiles(1).folder '/' myFiles(1).name];
                    NameSECFile = [curOutputSECDir '/' name '.sec'];
                end
                [x] = parseBinFile (NameFile, rP.dataType, rP.dataScalingIQ);   % parse the IQ data collected at Rx
                sig = x(1:rP.Ns);
                nG = 1;
                jj = (nG-1)*data.nFreqsperGrp+1:(nG)*data.nFreqsperGrp;
                [outputMat] = SEC_ALG(sig, data.D,  data.freq(data.indexFreqBins(jj)), false);
                [x2] = parseBinFile (NameSECFile, rP.dataType, 2^(-11));
                outputMat2     = buffer(x2, rP.Ns);
                norm(outputMat - outputMat2.')
            end
        else % Matlab processing
            [data, rP]= loadDataset( inputPathDir, data, rP);
            [rxSigN]= data.sig;
            if  rP.saveSigToMat
                save (NameSigMatFile, 'data', 'rP', '-v7.3');
            end

            % Run SEC + r.Window and sparseness of the received .mat matrix, save SEC.dat, SECSparse.dat
            if bitand( mode , 4 ) %  run SEC
                % Run SEC + r.Window and sparseness of the received .mat matrix
                [data] = processSEC(rxSigN, curOutputSECDir, curOutputSECSparseDir, NameSECMatFile, data, mode, rP);
            end
        end
    end
elseif bitand(mode , 64) % load SEC dataset

    % Create the SEC, SECSparse folders as well as sig.mat, SEC.mat files
    if ispc
        NameSECMatFile = [outputDir,  '\IntSEC.mat'];
    else
        NameSECMatFile = [outputDir,  '/IntSEC.mat'];
    end

    [data, rP]= loadDataset( inputPathDir, data, rP);
    % input Matlab sig, perform antenna integration for the given antsPerGroup
    [data, rP] = preprocess(data, [], rP, 3);

    if  rP.saveSigToMat
        save (NameSigMatFile, 'data', 'rP', '-v7.3');
    end
end

end

% Run SEC + r.Window and sparseness of the received .mat matrix
function [data] = processSEC(rxSigN, curOutputSECDir, curOutputSECSparseDir, NameSECMatFile, data, mode, rP)

[Ns, NmF]  = size(rxSigN);
rP.indSig = 1:Ns;
% input Matlab rx, run SEC, Max using reporting window (W) and sparsness (SP) at antenna level (params) - raturn matrices (SEC Full, Sparse
[data, rP] = preprocess(data, rxSigN, rP, 2);

if  rP.saveSigToMat
    save (NameSECMatFile, 'data', 'rP', '-v7.3');
end

if bitand( mode , 8 ) || bitand( mode , 16)
    for k = 1:NmF
        if ispc
            NameSECFile = [curOutputSECDir '\ch', sprintf('%02d',k), '_samples.dat'];
            NameSECSparseFile = [curOutputSECSparseDir '\ch', sprintf('%02d',k), '_samples.dat'];
        else
            NameSECFile = [curOutputSECDir '/ch', sprintf('%02d',k), '_samples.dat'];
            NameSECSparseFile = [curOutputSECSparseDir '/ch', sprintf('%02d',k), '_samples.dat'];
        end
        if rP.debugLog
            fprintf(' file name %s \n', NameSECFile);
        end

        if bitand( mode , 8 )
            % Full SEC output
            writeToFile (data.outputMat(:,k), NameSECFile, rP.dataType, rP.dataScalingSEC);
            % check if it saves correctly
            %[x2] = parseBinFile (NameSECFile, rP.dataType, rP.dataScalingSEC);
            % norm ( outputMat_(:) - x2)
            % Sparse SEC output
        end
        if bitand( mode , 16)
            writeToFile (data.STFT_sparse(:,k), NameSECSparseFile, rP.dataType, rP.dataScalingSEC);
            % check if it saves correctly
            %[x2Sprs] = parseBinFile (NameSECSparseFile, rP.dataType, rP.dataScalingSEC);
            % norm ( STFT_sparse_(:) - x2Sprs)
        end
    end
end
end




