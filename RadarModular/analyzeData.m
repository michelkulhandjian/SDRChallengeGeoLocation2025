% close all; clear all;

% Enable debug logs
rP.debugLog           = true;
isMoving              = false;
isFigSave             = true;
isVisual              = false;
isCompSNR             = false;
isMFMat               = false;
isEstimation          = false;
isDetection           = false;
isSignalLevel         = false;
isROC                 = false;
isCPP                 = false;

% Input directory
inputPathDir          = '/data/bin/';
movedPathDir          = '/data/bin/april18_2024';
outputPathDir         = '/home/michel.kulhandjian/radar/results/May_2024/SNR12';
logFile               = 'analyzeData.txt';
snrFile               = 'snrs.txt';
cppLogFile            = 'Cppresults.txt';
modes                 = [6]; % modes 1-signal visual, 2-MF output, 3- ROC plots, 4-compute SNR, 5-signal level 6-Estimation 7-Detection


clear folders;
% Dataset
%folders(1).name = 'bin';
% folders1(1).name = 'BIN1-A';
% folders1(2).name = 'BIN1-B';
% folders1(3).name = 'BIN2-A';
% folders1(4).name = 'BIN2-B';
% folders1(5).name = 'BIN3-A';
% folders1(6).name = 'BIN3-B';


[data, rP] = setParams('', inputPathDir, outputPathDir, rP);

%folders(1).isSignal = true;
%folders(2).isSignal = false;
h1MSs = []; h0MSs = [];
% ROC parameters
 rP.ROCsnr         = 25;
 rP.ROCnAntsperGrp = 48;
 rP.ROCWindow      = 30;
 rP.ROCNspelem     = 1;
 curWaveF       = 'BIN1-A';
 rP.ROCFig         = 1; % Plotting ROC figure type
 % 1. fix SNR, W, AntPerGrp, sp           - One ROC plot
 
 
% save it to res structure
res(1).SNR            = rP.ROCsnr;
res(1).ROCnAntsperGrp = rP.ROCnAntsperGrp;
res(1).ROCWindow      = rP.ROCWindow;
res(1).ROCNspelem     = rP.ROCNspelem;
res(1).ROCWaveform    = curWaveF;
res(1).ROCFig         = rP.ROCFig;
nameDetRes            = [outputPathDir '/det_res.mat'];
indSNR    = 1;
indAntInt = 1;
indW      = 1;

h1MSs = []; h0MSs = [];

if exist('movedPathDir','var')
    % fourth parameter does exist
    if isempty(movedPathDir)
        movedPathDir = [pwd];
    elseif ~exist(movedPathDir,'dir')
        mkdir(movedPathDir);
    end
else
    movedPathDir = [pwd];
end

% Open a results text file
fileID          = fopen([outputPathDir, '/', logFile] ,'a');
snrLogFile      = [outputPathDir, '/', snrFile];           

% Commands
%echo "password" | sudo -S ls -la /root/docker-test
sudoCmd            = 'echo "YOURPASSWORD" | sudo -S';
mkdir              = 'echo "YOURPASSWORD" | sudo -S -k mkdir -p ';
mv                 = 'echo "YOURPASSWORD" | sudo -S -k sudo mv ';
chown              = 'echo "YOURPASSWORD" | sudo -S -k chown -R michel.kulhandjian ';
%mkdir             = 'sudo mkdir -p ';
%mv                = 'sudo mv ';
%chown             = 'sudo chown -R michel.kulhandjian ';
%echo -e "YOURPASSWORD\n" | sudo -S yourcommand
pcreateenv         = 'python3 -m venv .venv';
psourceenv         = 'source .venv/bin/activate';
pinstallnumpy      = [psourceenv, ' ; ',  'pip install numpy'];
%pip install matplotlib
%pip install scipy
pplotter           = [psourceenv, ' ; ',  'python3 ../ploting_tools/data_parser_plotter_no_mtd.py -m -d '];

sourcePathDir = '/data/bin';
%rP.source = [sourcePathDir, '/code/radar-detector'];
%rP.build = [rP.source, '/build'];
rP.exper_file = [sourcePathDir, '/experiment2.json'];
rP.estimator  = [sourcePathDir, '/radar_estimator'];


NFFT      = 64;
NOVERLAP  = 60;
WINDOW    = 64;
C         = viridis(45);
numXticks = 4;


if isSignalLevel 
    if ~exist('.venv','dir')
        system(pcreateenv); system(psourceenv); system(pinstallnumpy);
    end
end

startA    = clock;
fprintf( 'Starting Analyzing dataset \n');
fprintf(fileID,' Starting Analyzing dataset \n' );

for indM = 1: length(modes)
    mode = modes(indM);
    clear folders;
    if mode == 1 % signal visual
        isMoving              = true;
        isVisual              = true;
        folders = folders1;
    elseif mode == 2 % MF outputs
        isMoving              = true;
        isVisual              = false;
        isMFMat               = true;
        folders(1).name = folders1(1).name;
    elseif mode == 3  %ROC plot
        isMoving              = false;
        isVisual              = false;
        isCompSNR             = false;
        isMFMat               = true;
        isROC                 = true;
        isSignalLevel         = false;
        folders(1).name = folders1(1).name;
        folders(2).name = '16-05-202418';  % Noise (ascii)
        folders(1).isSignal = true;
        folders(2).isSignal = false;
    elseif mode == 4  %Compute SNR
        isMoving              = false;
        isVisual              = true;
        isCompSNR             = true;
        isMFMat               = false;
        isROC                 = false;
        isSignalLevel         = false;
        folders(1).name = folders1(1).name;
        folders(2).name = '21-05-202416-25-36'; %'15-07-202413-59-03'; %02-07-202409-48-24'; %'21-05-202416-25-36'; % Noise (iq)
        folders(1).isSignal = true;
        folders(2).isSignal = false;
    elseif mode == 5  %Compute signal level
        isMoving              = false;
        isVisual              = false;
        isCompSNR             = false;
        isMFMat               = false;
        isROC                 = false;
        isSignalLevel         = true;
        folders(1).name = folders1(1).name;
    elseif mode == 6  %Estimate Radar parameters
        isMoving              = false;
        isVisual              = false;
        isCompSNR             = false;
        isMFMat               = false;
        isROC                 = false;
        isSignalLevel         = false;
        isEstimation          = true;
        folders               = folders1;
    elseif mode == 7  % Radar Detection
        isMoving              = false;
        isVisual              = false;
        isCompSNR             = false;
        isMFMat               = false;
        isROC                 = false;
        isSignalLevel         = false;
        isDetection           = true;
        folders               = folders1;
    end

    for dd = 1: length(folders)
        folderName   = folders(dd).name;
        pathDir      = [movedPathDir, '/', folderName];
        outputFolder = [outputPathDir, '/', folderName];

        if ~exist(outputFolder, 'dir')
            mkdirCmnd    = [mkdir, outputFolder];
            chownCmnd    = [chown, outputFolder];
            system(mkdirCmnd); system(chownCmnd);
        end

        if isMoving
            if ~exist(pathDir, 'dir')
                mkdirCommand = [mkdir, pathDir];
                chownCommand = [chown, pathDir];
                system(mkdirCommand); system(chownCommand);
            end
            mvCommand    = [mv, inputPathDir, '/*',  folderName, '*.dat ', pathDir, '/' ];
            mvascii      = [mv, inputPathDir, '/*',  folderName, '*.ascii ', pathDir, '/' ];
            %sudo mkdir -p 14-04-202408-12-18
            %sudo mv Det_14-04-202408-12-18*.dat 14-04-202408-12-18/
            system(mvCommand); system(mvascii);
        end

        if isVisual

            if rP.debugLog
                fprintf(' Visualizing IQ %s dataset \n', folderName);
            end
            fprintf(rP.fileID,' Visualizing IQ  %s dataset \n', folderName );

            [data, rP]   = generate_radar_dataset( pathDir, [], 32, rP);
            [rxsig]  = data.sig;
            [N, Nch] = size(rxsig);

            if isCompSNR
                if folders(dd).isSignal
                    sigP = sum(abs(rxsig(:)).^2)/numel(rxsig);
                else
                    noiseP = sum(abs(rxsig(:)).^2)/numel(rxsig);
                end
            else
                if rP.debugLog
                    fprintf(' Loading %s dataset \n', folderName);
                end
                fprintf(fileID,' Loading %s dataset \n', folderName);
                % Plot all the channels
                figSigAll = figure('visible','off');
                plot([1:N]/rP.fs*1e6, real(rxsig)); ylabel('Real Magnitude'); xlabel('Time (us)');
                %legendStrings = "Ch = " + string(rP.radioCh);
                %legend(legendStrings);
                %legend('Location', 'best');
                titleStrings = folderName + ", No Channels = " + string(Nch);
                title(titleStrings)
                if isFigSave
                    savefig(figSigAll,[outputFolder '/SignalAllChannels'  '.fig']);
                end
                saveas(figSigAll,[outputFolder '/SignalAllChannels'  '.png']);

                % select a channel
                indCh = randi([1 Nch],1,1);
                sig = rxsig(:,indCh);
                [~, indMax] = max(sig);
                indRange = indMax-3000:indMax+3000;
                if indRange(1) < 0
                    indRange = indRange+abs(indRange(1))+1;
                end
                if indRange(end) > N
                    indRange = indRange-abs(indRange(end))-1;
                end
                %[~, indR] = find (indRange > 0 & indRange<=N);
                sig = sig(indRange);
                figSig = figure('visible','off');
                plot([1:length(sig)]/rP.fs*1e6, real(sig)); ylabel('Real Magnitude'); xlabel('Time (us)');
                legendStrings = "Ch = " + string(indCh);
                legend(legendStrings);
                legend('Location', 'best');
                title(titleStrings);
                if isFigSave
                    savefig(figSig,[outputFolder '/SignalaChannel'  '.fig']);
                end
                saveas(figSig,[outputFolder '/SignalaChannel'  '.png']);

                % Plot Spectrogram
                figSpec = figure('visible','off');
                [~,Ftx,Ttx,PPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
                imagesc(Ttx/1e-3, Ftx/1e6, 10*log10(PPtx)); title("spectrogram");
                xticks(linspace(Ttx(1)/1e-3, Ttx(end)/1e-3, numXticks));
                xlabel('Time (ms)'); ylabel('Frequency (MHz)');
                set(gca,'YDir','normal'); colormap(C);
                title(titleStrings);
                if isFigSave
                    savefig(figSpec,[outputFolder '/spectrogram'  '.fig']);
                end
                saveas(figSpec,[outputFolder '/spectrogram'  '.png']);
            end
        end

        % Print SNRs
        if isSignalLevel
            if rP.debugLog
                fprintf(' Computing signal levels for  %s \n', folderName);
            end
            fprintf(rP.fileID,' Computing signal levels for  %s \n', folderName );
            plotterCmnd  = [pplotter, pathDir ' >> ', snrLogFile];
            %diary [outputPathDir, '/Snrs.txt']
            system(plotterCmnd);
            %diary off
        end

        if isMFMat
            if rP.debugLog
                fprintf(' Plotting MF outputs for %s \n', folderName);
            end
            fprintf(rP.fileID,' Plotting MF outputs for  %s \n', folderName );
            myFiles = dir(fullfile(pathDir,'*.ascii'));
            NmF = length(myFiles);

            if NmF > 0
                for k = 1:NmF
                    baseFileName = myFiles(k).name;
                    [ext, name] = fileparts(baseFileName);
                    if ispc
                        pathparts = strsplit(myFiles(k).folder,'\');
                        NameFile = [myFiles(k).folder '\' myFiles(k).name];
                    else
                        pathparts = strsplit(myFiles(k).folder,'/');
                        NameFile = [myFiles(k).folder '/' myFiles(k).name];
                    end
                    if rP.debugLog
                        fprintf(' Loading MF output %s \n', baseFileName);
                    end
                    fprintf(fileID,' Loading MF output %s \n', baseFileName );
                    matt = load(NameFile, '-ascii');

                    % ROC analysis
                    if isROC
                        if folders(dd).isSignal
                            h1MSs=[h1MSs matt.'];
                        else
                            h0MSs=[h0MSs matt.'];
                        end

                    else
                        figMFAll = figure('visible','off');
                        plot(matt); ylabel('MF values'); xlabel('number of MFs');
                        %legendStrings =  string(NameFile);
                        %legend(legendStrings);
                        %legend('Location', 'best');
                        %title(titleStrings);
                        title(extractBefore(baseFileName,".ascii") );
                        if isFigSave
                            savefig(figMFAll,[outputFolder 'DetectorMF' extractBefore(baseFileName,".ascii")  '.fig']);
                        end
                        saveas(figMFAll,[outputFolder 'DetectorMF' extractBefore(baseFileName,".ascii")  '.png']);
                    end
                end
            else
                if rP.debugLog
                    fprintf(' No MF output file is found in  %s \n', folderName);
                end
                fprintf(fileID,' No MF output file is found in  %s \n', folderName );
            end
        end

        if  isEstimation
            % estimation part
            if rP.debugLog
                fprintf(' Parameter estimating for  %s \n', folderName);
            end
            fprintf(rP.fileID,' Parameter estimating for  %s \n', folderName );
            if isCPP
                estLogFile = [outputFolder '/' cppLogFile ];
                %./radar_estimator --iq_sample_dir /data/data_share/May2024/30-05-202416-47-10 --json_config /data/bin/experiment.json
                cmdEst  = [rP.estimator, ' --iq_sample_dir ', pathDir, ' --json_config ',  rP.exper_file ' >> ', estLogFile];
                system(cmdEst);
            else
                runProcesses(pathDir, outputFolder, 'ESTIMATION', 4);
            end
        end
        if  isDetection
            % estimation part
            if rP.debugLog
                fprintf(' Radar Detection for  %s \n', folderName);
            end
            fprintf(rP.fileID,' Radar Detection for  %s \n', folderName );
            if isCPP
                %estLogFile = [outputFolder '/' cppLogFile ];
                %./radar_estimator --iq_sample_dir /data/data_share/May2024/30-05-202416-47-10 --json_config /data/bin/experiment.json
                %cmdEst  = [rP.estimator, ' --iq_sample_dir ', pathDir, ' --json_config ',  rP.exper_file ' >> ', estLogFile];
                %system(cmdEst);
            else
                runProcesses(pathDir, outputFolder, 'DETECTION', 4);
            end
        end
        if rP.debugLog
            fprintf(' finished looking into  %s \n', folderName);
        end
        fprintf(fileID,' finished looking into  %s \n', folderName);
    end

    if isCompSNR
        SNR = 10*log10(sigP/noiseP);
    end

    if isROC
        res(indSNR).sigpMSs{indAntInt}{indW} = h1MSs; %SNR, (AntPerGrp,W)
        res(indSNR).noisepMSs{indAntInt}{indW} = h0MSs;
        [res(indSNR).pd{indAntInt}{indW}, res(indSNR).pfa{indAntInt}{indW}, res.ROCthreshold] = roc_curve(h1MSs, h0MSs);
        plotROC ( res, rP )
        save (nameDetRes, 'res', 'rP');
        res.ROCthreshold
    end

end

time_s = (clock - startA)*[0 0 24*60^2 60.^[2 1 0]]';
elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
fprintf( 'Finished Analyzing Datasets %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
fprintf(fileID,'Finished Analyzing Datasets %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
fclose(fileID);