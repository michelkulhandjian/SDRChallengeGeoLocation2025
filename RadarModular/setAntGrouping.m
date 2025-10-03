function [data] = setAntGrouping( data)

% This will setup the antenna grouping parameters
% Given parameters are
% 1. Total number of active Antennas - this is dynamic depending on IQ
% 2. The number of antennas per group is set by the design
% 3. The number of frequency bins is set by the hardware constrain
% 4. The reporting window size W, set by the design
% Calculated parameters
% 1. The number of groups (nGroups)
% 2. The total number of frequency bins (totalFreqBins)
% 3. The indices (indexFreqBins)
% 4. The frequencies (freq)
% 5. The number of reporting windows in IQ stream

% Load fixed paramets
params

% Given parameters 
if isfield(data,'totalAntennas')
    totalAntennas = data.totalAntennas;
end

if isfield(data,'nAntsperGrp')
    nAntsperGrp = data.nAntsperGrp;
end

if isfield(data,'nFreqsperGrp')
    nFreqsperGrp = data.nFreqsperGrp;
end

if isfield(data,'WINDOW')
    WINDOW = data.WINDOW;
end


% compute antenna parameters
nGroups         = ceil(totalAntennas/nAntsperGrp);  
totalFreqBins   = nGroups*nFreqsperGrp;
indexFreqBins   = 1:1:totalFreqBins; %[1:2:totalFreqBins 2:2:totalFreqBins] ;
freq            = linspace(0,1, totalFreqBins).';  % Frequency bins
freq            = freq(1:totalFreqBins);

data.freq            = freq;
data.nGroups         = nGroups;
data.totalFreqBins   = totalFreqBins;
data.indexFreqBins   = indexFreqBins;
data.numW            = floor(Ns/WINDOW);










