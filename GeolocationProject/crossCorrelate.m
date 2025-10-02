function delay_es=crossCorrelate(s,sr,Nsample,fs)
    y=fftshift(ifft(fft(sr).*conj(fft(s))));        % Match Filter
    [~,I]=max(abs(y));                              % Find Maximum Index
    %delay_es=(I-Nsample/2-1)/fs;                    % Calculate Estimated Time Delay
    delay_es=(I-floor(Nsample/2+1))/fs;                    % Calculate Estimated Time Delay
    %figure; plot(abs(y))
end

