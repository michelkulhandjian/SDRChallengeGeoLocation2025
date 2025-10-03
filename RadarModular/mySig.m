function [N , sig] = mySig ( f0, fs, sweepBW, PW, PRF, numPulses, N)
    tau = 0;
    Ts = 1/fs;
    g1 = sweepBW/PW;
    T = fix(1./PRF*fs); % Pulse interval in samples
    %Np = fix(PW*fs);
    Np = ceil(PW*fs);
    nRange = [0:Np-1];
    sig = zeros(1, T);
    pulse= exp( 1i* ( 2*pi*f0*Ts*(nRange+tau) + pi*g1.*(Ts*(nRange+tau)).^2 ) );
    if N > Np
        pulse = [pulse, zeros(1, N-Np)];
        Np = N;
    end
    sig(1:Np) = pulse;
    sig = kron(ones(1,numPulses), sig);
    %sig = sig(1:Np).';
    sig = sig.';
end