% Preprocess Radar
function [data, rP] = preprocessRadar(data, sigA, rP, mode)

% =Preprocessing() has different mode as follows:
%  1. input Matlab rx, run SEC, Max using reporting window (W) and sparsness (SP) at antenna level - then perform antenna integration for the given antsPerGroup (params) 
%  2. input Matlab rx, run SEC, Max using reporting window (W) and sparsness (SP) at antenna level (params) - raturn matrices (SEC Full, Sparse)  
%  3. input Matlab SEC, perform antenna integration for the given antsPerGroup (params) 

if ~exist('mode','var')
    mode = 1;
end

if rP.debugLog
    %fprintf(' Start preprocessing process... \n');
end

% Write to results
%fprintf(rP.fileID,' Start preprocessing process... \n');

if mode == 2
    [Nsig, NmF] = size(sigA);
    data.totalAntennas   = NmF;
end

[data] = setAntGrouping( data);

numW        = data.numW ;
rWIN        = data.WINDOW;
nAnt        = data.nAntsperGrp;

data.numS = 0;
maxScheme = rP.maxScheme;  % 0-original, 1-antenna level, 2-antenna but with location

if (mode == 1 || mode == 2 || mode == 4) 
    sigA = sigA(rP.indSig, :);
    [Nsig, NmF] = size(sigA);
    rP.Ns  = size(sigA,1);
elseif ( mode == 3 )
    [Nsig, NmF] = size (data.STFT_sparse);
    rP.Ns  = Nsig/data.totalFreqBins;
end

if rWIN > Nsig
    rWIN = Nsig;
end
numW  = floor(Nsig/rWIN);
% This is to restrict larger number of processing
if 0
    NNN = 30;
    numW  = min(numW, NNN );
end
%--------
NNsig = fix(numW*rWIN);
if ( mode == 2 )
    data.outputMat = zeros (data.nFreqsperGrp*NNsig, NmF);
    data.STFT_sparse = zeros (data.nFreqsperGrp*NNsig, NmF);
elseif ( mode == 1 || mode == 3 || mode == 4 )
    data.STFT_sparse = zeros (data.totalFreqBins, NNsig);
    data.outputMat = zeros (data.totalFreqBins, NNsig);
end

for k = 1:NmF
    data.numS  = data.numS +1;
    if ( mode == 1 || mode ==2 || mode ==4) % Perform SEC + R.Windowing with delta-Max
        sig                    = sigA(1:NNsig,data.numS);
        data.data(:,data.numS) = sig;
    end 
    if ( mode == 1 || mode ==2 || mode == 4) % Perform SEC + R.Windowing with delta-Max
        if  rP.antInteg
            nG = ceil(data.numS/data.nAntsperGrp);
            jj = (nG-1)*data.nFreqsperGrp+1:(nG)*data.nFreqsperGrp;
            [outputMat] = SEC_ALG(sig, data.D,  data.freq(data.indexFreqBins(jj)), false);
        else
            % Output of SEC algorithm
            [outputMat] = SEC_ALG(sig, data.D, data.freq);
        end
        STFT_sparse = outputMat;
  %new      outputMatFull = outputMat;
        if maxScheme ==1 % 1-antenna level,
            STFT_sparse = [];
            for n = 1: numW
                [STFT_sparse]  = [STFT_sparse SparseP(double(outputMat(:,(n-1)*rWIN+1:n*rWIN)), data.num_sparse_elem, data.nFreqsperGrp ) ];
            end
    %new        outputMat = STFT_sparse;
        end
    end
    if ( mode == 3 )
        outputMat_  = data.sig(:,k);
        outputMat   = buffer(outputMat_, rP.Ns);
        outputMat   = outputMat.';
    end
    if ( mode == 1 || mode ==3 || mode == 4) % Antenna Integration
        if  rP.antInteg
            nG = ceil(data.numS/data.nAntsperGrp);
            index = (nG-1)*data.nFreqsperGrp+1:(nG)*data.nFreqsperGrp;
            curNumAntsperGrp = ((NmF-ceil(data.numS/nAnt)*nAnt)>=0)*nAnt+((NmF-ceil(data.numS/nAnt)*nAnt)<0).*(NmF-floor(data.numS/nAnt)*nAnt);
            %index = (data.numS-1)*data.nFreqsperGrp+1: data.numS*data.nFreqsperGrp;
            data.outputMat(index,:) = data.outputMat(index,:) + abs(outputMat)/curNumAntsperGrp;
            if ( mode == 1 && maxScheme ==1 ) % 1-antenna level
                data.STFT_sparse(index,:) = data.STFT_sparse(index,:) + abs(STFT_sparse)/curNumAntsperGrp;
            end
        else
            data.outputMat = data.outputMat + abs(outputMat)/NmF;
            if ( mode == 1 && maxScheme ==1 ) % 1-antenna level
                data.STFT_sparse = data.STFT_sparse + abs(STFT_sparse)/NmF;;
            end
        end
    elseif ( mode == 2)
        STFT_sparse_ = STFT_sparse.';
        outputMat_ = outputMat.';
        data.STFT_sparse(:,k) = STFT_sparse_(:);
        data.outputMat(:,k) = outputMat_(:);
    end
end

% Create Sparse choosing by choosing delta-Max
if ( mode == 1 || mode == 4)
    if maxScheme == 0   % 0-original,
        STFT_sparse = [];
        for n = 1: numW
            [STFT_sparse]  = [STFT_sparse SparseP(double(data.outputMat(:,(n-1)*rWIN+1:n*rWIN)), data.num_sparse_elem, data.totalFreqBins ) ];
        end
        %data.STFT_sparse = STFT_sparse(:, rP.indSig);
        data.STFT_sparse = STFT_sparse;
%new    else
 %new       data.STFT_sparse = data.outputMat;
    end
elseif ( mode == 3 )
    data.STFT_sparse = data.outputMat;
end



%fclose(rP.fileID);
end