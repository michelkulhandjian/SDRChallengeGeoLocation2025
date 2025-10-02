% Adding noise to received
function [y, noise] = addingAWGN ( x, reqSNR, mode)

if ~exist('mode','var')
    mode = 1;
end

%reqSNR = reqSNR + 17;
if (1)
%     threshold = 0.09 ; % percentage to take the signal for power calculations
%     indRxSig = abs(x(:,1)) > threshold*max(abs(x(:,1))) ;
%     %figure; plot(abs(x(:,1))); hold on; yline( threshold*max(abs(x(:,1))));
%     %figure; plot( real(x(indRxSig,1) ))
%     sig = x(indRxSig,1); %x(startInd:endInd);
    % 1e-5 for no channel applied
    indRxSig = abs(x(:)) > 1e-5 ;
    % 1e-3 for two-ray model & Rician (no channel as well)
    indRxSig = abs(x(:)) > 1e-3 ;
    if mode == 2
        indRxSig = abs(x(:)) > 25e-3 ;
        noise = x(~indRxSig);
        initNoisePower = sum(abs(noise(:)).^2)/numel(noise);
    end
    sig = x(indRxSig);
    sigPower = sum(abs(sig(:)).^2)/numel(sig); % linear
else
    fs              = 61.44e6; %30.72e6;
    [~, ~, ~, power] = obw(x(:), fs);
    sigPower = power/0.99;
    %sigPower = bandpower(x(:));
end
reqSNR = 10^(reqSNR/10);
noisePower = sigPower/reqSNR;
if mode == 2
    noisePower  = noisePower-initNoisePower;
end
noise = sqrt(noisePower/2)* (randn(size(x)) + 1i*randn(size(x)));
y = x + noise;
end