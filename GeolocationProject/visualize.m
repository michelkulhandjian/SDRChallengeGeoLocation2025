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

%% Plot Frequence response

% PSD 0
if 0
    figure; plot_psd(sig, rP.fs); xlabel('f (Hz)'); ylabel(' PSD'); grid on;
end

% PSD 1
if 1
    midfreq  = 0;
    window = blackman(NFFT);  % blackman window
    [pxx_1, f_vec] = pwelch(sig, window, [], NFFT, rP.fs, 'centered', 'power');
    % Convert power to dB scale
    pxx_dB = 10*log10(pxx_1);

    % Adjust frequency axis for center frequency
    f_axis = f_vec + midfreq;
    figure;
    plot(f_axis/1e6, pxx_dB);  % Convert frequency to MHz for better visualization
    grid on;
    axis tight
    xlabel('Frequency (MHz)');
    ylabel('Power (dB)');
    title ('Spectrum')
    %title(sprintf('Spectrum of Chunk %d', n));  % Dynamic title with chunk number
    % savefig(sprintf('%s/Spectrum_Chunk_%d.fig', iq_folder, n));  % Save as .fig
    % saveas(gcf, sprintf('%s/Spectrum_Chunk_%d.png', iq_folder, n));  % Save as .png
    % Set figure and axes properties to ensure consistency
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 600, 400]);  % Fixed figure size
    set(gca, 'LooseInset', get(gca, 'TightInset'));  % Minimize padding
end

if 0
    fs = rP.fs;
    f = [-fs/2:fs/NFFT:fs/2-.0001];
    iq = sig;
end

% PSD 2
if 0
    px=pwelch(iq,blackman(NFFT),32, NFFT, rP.fs, 'twosided');
    figure;plot(f,10*log10(abs(fftshift(px))));grid on;
    title('Spectrum 1')
    %title(strcat('Spectrum of - ', name));
    xlabel ('MHz'); ylabel ('dB');
end

% PSD 3
if 0
    P_size = NFFT;
    px=pwelch(iq,blackman(P_size),0, P_size, rP.fs, 'twosided');
    f = [-fs/2:fs/P_size:fs/2-.1*fs/P_size];
    figure;
    plot(f,10*log10(abs(fftshift(px))));
    grid on;
    title('Spectrum 2')
    %title(strcat("Spectrum of => ","" ,name));
    xlabel ('MHz'); ylabel ('dB');
end