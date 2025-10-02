% Parsing the IQ data
function [x] = parseBinFile (filename, dataType, dataScaling)
fid = fopen(filename, 'r');
data = fread(fid, inf, dataType);
fclose(fid);
I = data(1:2:end);
Q = data(2:2:end);
num_samps = length(I);
I = I(1:num_samps)*dataScaling;
Q = Q(1:num_samps)*dataScaling;
x = complex(I, Q);

end

