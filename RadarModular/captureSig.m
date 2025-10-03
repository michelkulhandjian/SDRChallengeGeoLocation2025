function [fig1, fig2, fig3] = captureSig (sig, rP, data)

C       = viridis(45);
numXticks = 4;
if rP.Spec
    outputMat = sig;
    fig1 = figure('visible','off');
    imagesc(rP.timeVec/1e-3, data.freq*rP.fs/1e6, abs(double(outputMat)));
    xlabel('Time (ms)'); ylabel('Frequency (Mhz)'); title('SEC output');
    set(gca,'YDir','normal'); colormap(C);

    fig2 = 0;
    fig3 = 0;

else
    fig1 = figure('visible','off');
    [outputMat] = SEC_ALG(sig, data.D, data.freq);
    imagesc(rP.timeVec/1e-3, data.freq*rP.fs/1e6, abs(double(outputMat)));
    xlabel('Time (ms)'); ylabel('Frequency (Mhz)'); title('SEC output');
    set(gca,'YDir','normal'); colormap(C);

    fig2 = figure('visible','off');
    plot (real(sig));
    xlabel('Samples'); ylabel('Real Mag.'); title('SIgnal');

    fig3 = figure('visible','off');
    % Plot Radar Pulse
    NFFT = 64;
    NOVERLAP = 60;
    WINDOW = 64;
    [~,Ftx,Ttx,PPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
    imagesc(Ttx/1e-3, Ftx/1e6, 10*log10(PPtx)); title("spectrogram");
    xticks(linspace(Ttx(1)/1e-3, Ttx(end)/1e-3, numXticks));
    xlabel('Time (ms)'); ylabel('Frequency (MHz)');
    set(gca,'YDir','normal'); colormap(C);
end
end