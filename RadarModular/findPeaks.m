% Finding peaks
function [cpeaks, TEst, PWEst] = findPeaks ( sig)

PLOT_FIG = true;
margin = 60;
adj    = 3.4372; %1.08;
[cpeaks, ~] = find (abs(real(sig)) > mean(abs(real(sig)))*adj);

consInd = diff(cpeaks)==1;
consM = find([false;consInd]~=[consInd;false]);
consS   = consM(2:2:end)-consM(1:2:end-1);
sInd = find(consS>=margin);
if ~isempty(sInd)
    indP = floor( (cpeaks(consM(2*sInd-1)) + cpeaks(consM(2*sInd)) )/2);
    PWEst = mean( cpeaks(consM(2*sInd)) - cpeaks(consM(2*sInd-1)));

    if ~isempty(indP) && length(indP) > 1
        TEst = mean(diff(indP));
    else
        TEst = length(sig);
    end
    cpeaks = indP';
    if PLOT_FIG
        figure; plot(sig)
        hold on; xline(indP);
    end
else
    TEst = 600;
    PWEst = 600;
    cpeaks = [1:600];
end
end

