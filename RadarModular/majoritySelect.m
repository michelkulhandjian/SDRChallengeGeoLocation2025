function [inxAnomaly, indM, EstM] = majoritySelect (X, perc)

if ~exist('perc','var')
    perc = 0.10;  % percentage removing from the X
end

inxAnomaly  = [];
%[N,edges]   = histcounts(X);
[N,edges]   = histcounts(X, 100);
posInx      = find(N);
indM        = 1;
if length( unique(N(posInx)) ) ~= 1
    [~, maxInx] = maxk(N,1);
    [~, inx2 ] = find (X >= edges(maxInx) & X <= edges(maxInx+1));
    maxX       = max(X(inx2));
    EstM       = mean(X(inx2));
    indM       = inx2(1);
    if 0
        frac = N/sum(N);
        [~, inx ]  = find (frac < perc);
        inxBins = intersect(posInx, inx);
        for ii = 1: length(inxBins)
            [~, inx2 ] = find (X >= edges(inxBins(ii)) & X <= edges(inxBins(ii)+1));
            inxAnomaly  = [inxAnomaly inx2];
        end
    else
        [~,inxAnomaly] = find ((X<maxX*(1-perc)) | (X>maxX*(1+perc)));
    end
else
    EstM = mean(X);
end

