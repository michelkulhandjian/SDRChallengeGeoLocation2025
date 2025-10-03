% estimate Waveform
function [out] = selectWaveform(pw, pri, bw, isSelect, rP)

distPW = ([1 2 50 100 10 50 50 10]*1e-6 - pw)/1e-6;
%distPri = ([50 40 10000 6250 20000 10000 10000 20000]*1e-6 - pri)/1e-6;
% New one
distBW = ([1 8 20 30 5 15 20 5]*1e6 - bw)/1e6;
distPri = ([50 33.3 5000 5000 10000 6666.7 10000 20000]*1e-6 - pri)/1e-6;

if isSelect
    dist = (2*abs(distPW) + abs(distPri))/2;
else
    dist = (abs(distPW) + abs(distBW) + abs(distPri));
end
[~, ind] = min(dist);
ALL = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B', 'OTHER-A', 'OTHER-B'};
mode = ALL{ind(1)};
if isSelect
    rP.mode =mode;
    rP = setRadarRange( rP);
    out = rP;
else
    out = mode;
end

end