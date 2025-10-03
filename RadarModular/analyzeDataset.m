% close all; clear all;

% Enable debug logs
rP.debugLog           = true;
isMoving              = true;
isFigSave             = false;
isVisual              = true;
isMFMat               = false;
isEstimation          = false;
isSignalLevel         = false;

% Input directory
inputPathDir          = '/data/bin/';
%movedPathDir          = '/dumps/IQDumpFieldTest/tests_iq_dumps/april17_2024/';
movedPathDir          = '/data/bin/april18_2024';
movedPathDir          = '/data/data_share/RadarGoldenDataset2';
outputPathDir         = '/data/bin/april18_2024';
outputPathDir         = '/home/michel.kulhandjian/radar/results/april23_2024';
logFile               = 'analyzeData.txt';
snrFile               = 'snrs.txt';


% Dataset
%folders(1).name = 'bin';
%folders(1).name = '18-04-202409-10-19';
% folders(2).name = '18-04-202414-54-57';
% folders(3).name = '18-04-202414-53-17';
% folders(4).name = '18-04-202414-53-04';
% folders(5).name = '18-04-202414-51-11';
% folders(6).name = '18-04-202414-50-58';
% folders(7).name = '18-04-202414-49-38';
% folders(8).name = '18-04-202414-49-08';
% folders(9).name = '18-04-202414-48-55';
% folders(10).name = '18-04-202413-02-50'; %out
% folders(11).name = '18-04-202414-48-43';

[data, rP] = setParams('', inputPathDir, outputPathDir, rP);

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
chown              = 'echo "YOURPASSWORD" | sudo -S -k chown -R YOURUSER ';
%mkdir             = 'sudo mkdir -p ';
%mv                = 'sudo mv ';
%chown             = 'sudo chown -R michel.kulhandjian ';
%echo -e "YOURPASSWORD\n" | sudo -S yourcommand
pcreateenv         = 'python3 -m venv .venv';
psourceenv         = 'source .venv/bin/activate';
pinstallnumpy      = [psourceenv, ' ; ',  'pip install numpy'];
%pip install matplotlib
%pip install scipy
pplotter           = [psourceenv, ' ; ',  'python3 ../ploting_tools/data_parser_plotter.py -m -d '];


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

    % Print SNRs
    if isSignalLevel
        if rP.debugLog
            fprintf(' Computing SNRs for  %s \n', folderName);
        end
        fprintf(rP.fileID,' Computing SNRs for  %s \n', folderName );
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

                figMFAll = figure('visible','off');
                plot(matt); ylabel('MF values'); xlabel('number of MFs');
                %legendStrings =  string(NameFile);
                %legend(legendStrings);
                %legend('Location', 'best');
                %title(titleStrings);
                title(extractBefore(baseFileName,".ascii") );
                if isFigSave
                    savefig(figMFAll,[outputFolder '/DetectorMF' extractBefore(baseFileName,".ascii")  '.fig']);
                end
                saveas(figMFAll,[outputFolder '/DetectorMF' extractBefore(baseFileName,".ascii")  '.png']);
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
        runProcesses(pathDir, outputFolder, 'ESTIMATION', 4);

    end
    if rP.debugLog
        fprintf(' finished looking into  %s \n', folderName);
    end
    fprintf(fileID,' finished looking into  %s \n', folderName);
end



time_s = (clock - startA)*[0 0 24*60^2 60.^[2 1 0]]';
elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
fprintf( 'Finished Analyzing Datasets %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
fprintf(fileID,'Finished Analyzing Datasets %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
fclose(fileID);