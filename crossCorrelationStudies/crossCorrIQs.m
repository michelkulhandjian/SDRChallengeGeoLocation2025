%% Cross correlation of IQs

dataType        = 'int16';  % Set the data read type (e.g., 'int16', 'float32', etc.)
dataScalingIQ   = 2^(-15);  % This is scaling for IQ data
offset1         = 0;        % offset for rx1
offset2         = 0;        % offset for rx2
Nbuffer         = 1000;     % Buffer or IQ to be processed
rP.fs           = 61.44e6;  % Rx sampling rate;
c               = physconst('LightSpeed');

NameFile1 = 'receiver1.dat';
[x1]  = parseBinFile (NameFile1, dataType, dataScalingIQ); 

NameFile2 = 'receiver21.dat';
[x2]  = parseBinFile (NameFile2, dataType, dataScalingIQ); 

%figure; plot(real(x1(1:10000))); hold on;
%plot(real(x2(1:10000)));

rx1 = x1([1:Nbuffer]+offset1).';
rx2 = x2([1:Nbuffer]+offset2).';

visualize (rx1, rP, 'rx1 time domain', 'rx1');
visualize (rx2, rP, 'rx2 time domain', 'rx2');

x = linspace(0,1,length(rx1));

[corr0, delaySample0] = crossCorrelationAll (rx1, rx2, 0);
[corr1, delaySample1] = crossCorrelationAll (rx1, rx2, 1);
[corr2, delaySample2] = crossCorrelationAll (rx1, rx2, 2);
[corr3, delaySample3] = crossCorrelationAll (rx1, rx2, 3);
[corr4, delaySample4] = crossCorrelationAll (rx1, rx2, 4);
[corr5, delaySample5] = crossCorrelationAll (rx1, rx2, 5);
[corr6, delaySample6] = crossCorrelationAll (rx1, rx2, 6);

% Print results
fprintf("max corr0: %d, corr1: %d, corr2: %d, corr3: %d, corr4: %d, corr5: %d, corr6: %d\n", delaySample0, delaySample1, delaySample2, delaySample3, delaySample4, delaySample5, delaySample6);



%% Plot the results

% ---------- Plotting ----------
figure;

subplot(2,1,1);
plot(x, real(rx1), 'b'); hold on;
plot(x, real(rx2), 'r');
grid on; xlim([0 1]); legend(  'rx1', 'rx2' );
title("Time Domain: rx1 and rx2");

Ncorr = length(corr0);

subplot(2,1,2);
plot(0:length(corr0)-1, real(corr0), '.-'); hold on;
plot((0:length(corr1)-1), real(corr1), '.-');
plot((0:length(corr2)-1), corr2, '--');
plot((0:length(corr4)-1), corr4, '-*');
plot((0:length(corr5)-1), corr5, '-x');
grid on;   legend(  'FFT-based', 'xcorr', 'FFT-based 1',    'FFT-based 2',  'GCC-PHAT' );
xlim([0 Ncorr]);
title("Different cross-correlation");

figure; plot(0:length(corr3)-1, corr3, '-*'); hold on;
plot((0:length(corr6)-1), real(corr6), '.-');
grid on;   legend(  'FFT-based zero padding',   'FFT-based zero padding2' );
title("Different cross-correlation");

% need to find out what is the best cross correlation function we can use?
tdoa21 = delaySample5/rP.fs;
rangeDiff21 = tdoa21*c;   % in meters - does this distance make sense?
