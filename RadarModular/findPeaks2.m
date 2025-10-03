% FindPeaks2
% Finding peaks
function [rpeaks, peaks, TEst, TEst1, PWEst] = findPeaks2 (sig, rP)
%Np = rP.fs
margin              = 170;
pulseOffset =0;

%[cpeaks, ~] = find ((real(sig)) > max(real(sig))*rP.percent_peak);
[cpeaks, ~] = find (abs(real(sig)) > max(abs(real(sig)))*rP.percent_peak);
%corr2 = conv(conj(flipud(sig(1:Np))), sig);
%cpeaks = find(abs(corr2) > percent_peak*max(abs(corr2)));
if ~isempty(cpeaks)
    [ind, ~] = find ( [0;diff(cpeaks)] > margin) ;
    %xline(cpeaks(ind(1:10))

    Lcpeaks = length(cpeaks);
    L = length(sig);
    startInd = [1;ind];   % Please check here, do we need [1;ind] only?
    endInd = [ind-1;Lcpeaks];

    [indm, ~] =find (startInd >Lcpeaks);
    if ~isempty(indm)
        startInd(indm) = Lcpeaks;
    end
    [indm, ~] =find (endInd <1);
    if ~isempty(indm)
        endInd(indm) = 1;
    end
    peaks = [cpeaks(startInd) cpeaks(endInd)];
    PW = [];
    T1 = [];

    curPulses = length(startInd);
    if curPulses > (rP.numPulses+pulseOffset)
        curPulses = rP.numPulses +pulseOffset;
    end
else
    curPulses = 0;
    peaks    = 0;
end
rpeaks = 1;
if curPulses > 1+pulseOffset
    for nn = 1+pulseOffset:curPulses
        % Compute PW
        curPW = abs(cpeaks(endInd(nn)) - cpeaks(startInd(nn)));
        PW = [PW curPW];

        % Compute T
        sumT = cpeaks(startInd(nn))+cpeaks(endInd(nn));
        rpeaks(nn) = fix(sumT/2);
        % Compute T1
        if nn ~= 1+pulseOffset
            T1 = [T1  (cpeaks(startInd(nn)) - cpeaks(startInd(nn-1))) ];
        end
    end
    % Estimating T
    TEst = diff(rpeaks);
    [TEst, indMM] = selectMajoritySimilarValues (TEst, margin);
    rpeaks = rpeaks(indMM);

    % Estimating T1
    %[~,indMM1] = find ( abs(mode(T1) - T1) < margin);
    %if ~isempty(indMM1)
    %TEst1 = mean(T1(indMM1));
    %else
    %TEst1 = mean(T1);
    %end
    [TEst1, indMM1] = selectMajoritySimilarValues (T1, margin);

    % Estimating PW
    PWEst = max(PW);
    %PWEst = mean(PW);
else
    TEst = length(sig);
    TEst1 = TEst;
    PWEst = TEst;
end
end

function [TEst1, indMM1] = selectMajoritySimilarValues ( T1, margin)

indMM1 = 1:length(T1);
TEst1  = mean(T1);
compT1 = T1 - mean(T1);
[~,ind] = find ( abs(compT1) < margin);
if (isempty(ind) ) || (length(ind) ~= length(indMM1) )

    [~,indMPos] = find ( compT1 >0 );
    [~,indMNeg] = find ( compT1 <0 );
    [~,indMZero] = find ( compT1 == 0 );
    if ~isempty(indMPos)
        if ~isempty(indMZero)
            if ~isempty(indMNeg) && abs(mean(compT1(indMPos))) <= abs(mean(compT1(indMNeg)))
                indMPos = [indMPos indMZero];
            else
                indMNeg = [indMNeg indMZero];
            end
        end

        if ~isempty(indMNeg)  && length(indMNeg) > length(indMPos)
            indMM1 = [indMNeg];
        else
            indMM1 = [indMPos];
        end
        TEst1 = mean(T1(indMM1));
    elseif ~isempty(indMNeg)
        if ~isempty(indMZero)
            if ~isempty(indMPos) && abs(mean(compT1(indMPos))) <= abs(mean(compT1(indMNeg)))
                indMPos = [indMPos indMZero];
            else
                indMNeg = [indMNeg indMZero];
            end
        end

        if ~isempty(indMPos)  && length(indMNeg) > length(indMPos)
            indMM1 = [indMNeg];
        else
            indMM1 = [indMPos];
        end
        TEst1 = mean(T1(indMM1));
    end
    indMM1 = unique([indMM1 indMM1+1]);
    indMM1 = indMM1(2:end-1);
end
end