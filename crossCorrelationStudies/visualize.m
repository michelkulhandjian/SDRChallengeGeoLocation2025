function visualize (sig, rP, legendStrings, titleStrings)
% plot signal
%sig = sig(indRange);
%figSig = figure('visible','off');
figure;
plot([1:length(sig)]/rP.fs*1e6, real(sig)); ylabel('Real Magnitude'); xlabel('Time (us)');
%legendStrings = "Ch = " + string(indCh);
legend(legendStrings);
legend('Location', 'best');
title(titleStrings);

NFFT      = 64;
NOVERLAP  = 60;
WINDOW    = 64;
C         = viridis(45);
numXticks = 4;
% Plot Spectrogram
%figSpec = figure('visible','off');
figure;
[~,Ftx,Ttx,PPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
imagesc(Ttx/1e-3, Ftx/1e6, 10*log10(PPtx)); title("spectrogram");
xticks(linspace(Ttx(1)/1e-3, Ttx(end)/1e-3, numXticks));
xlabel('Time (ms)'); ylabel('Frequency (MHz)');
set(gca,'YDir','normal'); colormap(C);
title(titleStrings);
