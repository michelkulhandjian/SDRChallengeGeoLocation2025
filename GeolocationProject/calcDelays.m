function [tdoaEst] = calcDelays(sig, rP)


[L, N]     = size(sig);
if rP.analyzeCross
    Nbuffer    = 1000;     % Buffer or IQ to be processed
    offset1    = 0;        % offset for rx1
    offset2    = 0;        % offset for rx2
    if Nbuffer > L
        Nbuffer = L;
    end
    rx1 = sig([1:Nbuffer]+offset1,1).';
    visualize (rx1, rP, 'rx1 time domain', 'rx1');
    x = linspace(0,1,length(rx1));

    [tDelay]=calc_Delay(rP.bPos, rP.pPos);
    tdoaData = tDelay(2:end)-tDelay(1);
end



for n = 2: N
    tdoaEst(n-1)=crossCorrelate(sig(:,1),sig(:,n),L,rP.fs);
    
    %[corr, delaySample] = crossCorrelationAll (sig(:,1),sig(:,n), rP.crossMode);
    %tdoaEst(n-1)= delaySample/rP.fs;

    if rP.analyzeCross
        rx2 = sig([1:Nbuffer]+offset2,n).';

        [corr0, delaySample0] = crossCorrelationAll (rx1,rx2, 0);
        [corr1, delaySample1] = crossCorrelationAll (rx1,rx2, 1);
        [corr2, delaySample2] = crossCorrelationAll (rx1,rx2, 2);
        [corr3, delaySample3] = crossCorrelationAll (rx1,rx2, 3);
        [corr4, delaySample4] = crossCorrelationAll (rx1,rx2, 4);
        [corr5, delaySample5] = crossCorrelationAll (rx1,rx2, 5);
        [corr6, delaySample6] = crossCorrelationAll (rx1,rx2, 6);

        Ncorr = length(corr0);

        % Print results
        fprintf("max corr0: %d, corr1: %d, corr2: %d, corr3: %d, corr4: %d, corr5: %d, corr6: %d\n", delaySample0, delaySample1, delaySample2, delaySample3, delaySample4, delaySample5, delaySample6);
        fprintf("Delay corr0: %d, corr1: %d, corr2: %d, corr3: %d, corr4: %d, corr5: %d, corr6: %d\n", delaySample0/rP.fs, delaySample1/rP.fs, delaySample2/rP.fs, delaySample3/rP.fs, delaySample4/rP.fs, delaySample5/rP.fs, delaySample6/rP.fs);
        
        fprintf("True Delay: %d \n", tdoaData(n-1) );

        
        visualize (rx2, rP, ['rx', num2str(n), ' time domain'], ['rx',num2str(n)]);

        figure;

        subplot(2,1,1);
        plot(x, real(rx1), 'b'); hold on;
        plot(x, real(rx2), 'r');
        grid on; xlim([0 1]); legend(  'rx1', ['rx', num2str(n)]);
        tt = sprintf('Time Domain: rx1 and rx%d', n);
        %title("Time Domain: rx1 and rx2");
        title(tt);

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

    end

end
