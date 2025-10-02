% loads the dataset for a given inputPath then converts to MAT file
function [data, rP]= loadDataset( inputPathDir,data, rP)

if rP.debugLog
    fprintf(' Start data loading process... \n');
end
% Write to results
fprintf(rP.fileID,' Start data loading process... \n');

if ~rP.isSingleFile
    myFiles = dir(fullfile(inputPathDir,'**','*.dat')); %gets all mat files in struct
    NmF = length(myFiles);
else
    NmF = 1;
end

% the maximum length of signal per file
Ns       = ceil(rP.Ts*rP.fs);
%data.sig  = zeros(Ns,1);
%clear data.sig ;
data = rmfield(data,'sig');

if NmF > 0
    data.totalAntennas   = NmF;
    for k = 1:NmF
        if ~rP.isSingleFile
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
                fprintf(' file name %s \n', baseFileName);
            end
        end
		chNum = regexp(baseFileName,'\d*','Match');
        [x]  = parseBinFile (NameFile, rP.dataType, rP.dataScalingIQ);   % parse the IQ data collected at Rx
		
		rP.channels = [rP.channels str2num(chNum{1})];  % Signal Index
        if numel(x) > Ns
            sig  = x(1:Ns);
        else
            sig  = x;
        end
        if sum(x) ~= 0
		    data.numS     = data.numS +1;
            data.sig(:, data.numS) = sig;
        end
    end

    %save (rP.nameProcSECData, 'data', 'rP', '-v7.3');
    if rP.debugLog
        fprintf(' Completed Data loading Process... \n');
    end

    % Write to results
    fprintf(rP.fileID,' Completed Data loading Process... \n');

else
    if rP.debugLog
        fprintf(' No Files are found... \n');
    end
    fprintf(rP.fileID,' No Files are found...  \n' );
end
%fclose(rP.fileID);
end