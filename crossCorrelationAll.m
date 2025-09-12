function [corr, delaySample] = crossCorrelationAll (tx, rx, method_flag)


if method_flag == 0
    %%% FFT based cross correlation
    Ntx = length(tx);
    corr=fftshift(ifft(fft(tx).*conj(fft(rx))));        % Match Filter
    [~,I]=max(abs(corr));                              % Find Maximum Index
    delaySample =(I-floor(Ntx/2+1));                    % Calculate Estimated Time Delay

elseif method_flag == 1
    % ---------- Cross correlation using convolution (xcorr) ----------
    % Python's signal.correlate(mode="same") ≈ centered version of xcorr
    corr_full = xcorr(tx, rx);
    mid = ceil(length(corr_full)/2);
    Ntx = length(tx); N1 = floor(length(tx)/2);
    corr = corr_full(mid - N1 : mid + floor((Ntx-1)/2));
    [~, idx_max] = max(corr);
    delaySample = idx_max - N1;

elseif method_flag == 2
    fft_decim = 1;        % Sample drop ratio before doing FFT (1 -> nothing dropped)
    N1 = floor(length(tx)/2);
    % ---------- FFT/IFFT based cross correlation ----------
    TX = fft(tx(1:fft_decim:end));
    RX = fft(rx(1:fft_decim:end));
    CORR = TX .* conj(RX);
    corr = fftshift(real(ifft(CORR)));
    [~, idx_max_fft] = max(corr);
    delaySample = (idx_max_fft - N1) * fft_decim;

elseif method_flag == 3
    fft_decim = 1;        % Sample drop ratio before doing FFT (1 -> nothing dropped)
     N1 = floor(length(tx));
    TX1 = fft([tx(1:fft_decim:end) zeros(size(tx))]);
    RX1 = fft([rx(1:fft_decim:end) zeros(size(rx))]);
    CORR1 = TX1 .* conj(RX1);
    corr = fftshift(real(ifft(CORR1)));
    [~, idx_max_fft1] = max(corr);
    delaySample = (idx_max_fft1 - N1) * fft_decim;

elseif method_flag == 4
    Ntx = length(tx);
    p = nextpow2(Ntx);
    FFTLength=2^p;
    N1 = floor(Ntx/2);
    %Rxx = xcorr(y1);
    %Ryy = xcorr(y2);
    Rxy = xcorr(tx,rx);
    %Sxx = fft(Rxx,FFTLength);
    %Syy = fft(Ryy,FFTLength);
    Sxy = fft(Rxy,FFTLength);
    % Unfiltered Correlation (plain correlation)
    W = ones(1,FFTLength);
    % Apply the filter
    R = Sxy.*W;
    % Obtain the GCC
    corr = fftshift(real(ifft(R)));
     [~, idx_max] = max(corr);
    delaySample = idx_max - N1;

elseif method_flag == 5

    Ntx = length(tx);
    p = nextpow2(Ntx);
    FFTLength=2^p;
    N1 = floor(Ntx/2);
    Rxy = xcorr(tx,rx);
    Sxy = fft(Rxy,FFTLength);
    % PHAT Filter
    W = 1./abs(Sxy);
    % Apply the filter
    R = Sxy.*W;
    % Obtain the GCC
    corr = fftshift(real(ifft(R)));
    [~, idx_max] = max(abs(corr));
    delaySample = idx_max - N1;

elseif method_flag == 6

    len_tx = length(tx);
    len_rx = length(rx);
    N1 = floor(length(tx));
    N_padded = 2^nextpow2(len_tx + len_rx- 1);
    tx_padded = [tx, zeros(1, N_padded - len_tx)];
    rx_padded = [rx, zeros(1, N_padded - len_rx)];

    fft_tx = fft(tx_padded);
    fft_rx = fft(rx_padded);
    R_xy_fft = fft_tx .* conj(fft_rx);
    cross_correlation = fftshift(ifft(R_xy_fft));

    if 0
        corr = cross_correlation;
    else
        % The full length of the linear cross-correlation is len_x + len_y - 1
        linear_corr_length = len_tx + len_rx - 1;
        corr = cross_correlation(1:linear_corr_length);
    end

    [~, idx_max] = max(corr);
    delaySample = idx_max - N1;

    % % Display the result
    % disp('Cross-correlation with zero-padding:');
    % disp(real(cross_correlation));
    %
    % % Compare with MATLAB's built-in xcorr function (which also uses zero-padding internally)
    % disp('Cross-correlation using xcorr:');
    % disp(xcorr(x,y));

else
    disp('This option does not exist');

end

end


function xcor1(inputSize, iq1, iq2, scale, lag)
    % Initialize variables
    xcor = zeros(1, 2 * lag + 1); % Preallocate xcor array
    for m = -lag:lag
        corr_real = 0; 
        corr_img = 0;
        end1 = (m < 0) * (inputSize + m) + (m >= 0) * (inputSize - m);
        for n = 0:end1-1
            indx = (m < 0) * n + (m >= 0) * (n + m);
            indy = (m < 0) * (n - m) + (m >= 0) * n;

            real1 = scale * iq1(2 * indx + 1);
            img1 = scale * iq1(2 * indx + 2);
            real2 = scale * iq2(2 * indy + 1);
            img2 = scale * iq2(2 * indy + 2);

            corr_real = corr_real + (real1 * real2 + img1 * img2);
            corr_img = corr_img + (real1 * img2 - real2 * img1);
        end

        xcor(m + lag + 1) = 10.0 * log10(sqrt(corr_real^2 + corr_img^2) / inputSize);
    end
end