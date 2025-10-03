function [outputMat] = SEC_ALG(rxSig, D, freq, isFixedPoint)

if nargin == 3   % if the number of inputs equals 3
    isFixedPoint = true; % then make the isFixedPoint value, equal to  default value, true.
end

if isFixedPoint % Fixed point implementation
    %% Sine Wave (detection)
    nWord          = 16;   % Word length for fixed-point implementation (nWord - nFract = Int. length)
    nFract         = 11;   % fractional lenght, can be modified but output mat files should be changed as well!
    sidx           = 0;    % Index
    numFreqs       = length(freq);    % length of the frequency bins
    N              = length(rxSig);   % length of the signal
    nn             = 0:1:N-1;         % row vector for n
    %   fi(Data, Signed, WordLength, FractionLength)
    rxSig          = fi(rxSig, 1, nWord, nFract);

    % Initialization
    sineWave    = fi ( double( exp(-1i*2*pi*freq*nn) ), 1, nWord, nFract);
    in_sine     = sineWave(:, 1:N);  % sine "buffer"
    outputMat   = fi( zeros(numFreqs, N), 1, nWord, nFract);
    delayBuffer = fi( zeros(numFreqs, D), 1, nWord, nFract);
    sineBuffer  = fi( zeros(numFreqs,D), 1, nWord, nFract);   % For simplicity (delay buffer for sine)
    outVal      = fi( 0, 1, nWord, nFract);
    outVec      = fi( [], 1, nWord, nFract);
    while ( sidx < N)
        sidx = sidx + 1;
        % New samples
        newSamp    = fi( in_sine(:,sidx) * rxSig(sidx), 1, nWord, nFract);
        buffPtr    = mod(sidx-1, D) + 1;
        buffer_out = fi( delayBuffer(:,buffPtr) .* sineBuffer(:,buffPtr), 1, nWord, nFract) ; % previous iter
        sineBuffer(:,buffPtr)  = in_sine(:,sidx);
        delayBuffer(:,buffPtr) = rxSig(sidx) ;
        % Out
        outVal     = fi ( outVal + newSamp - buffer_out, 1, nWord, nFract);
        %  These last two lines might be uneccessary, we can come back to this
        %  at the later stage
        %buffPtrAdj = mod(sidx, D) + 1;
        %outV       = fi( outVal.*conj(sineBuffer(:,buffPtrAdj)), 1, nWord, nFract);
        outV       = outVal;
        outVec     = [outVec outV];
    end
    outputMat = outVec;
else % Floating point implementation    
    sidx           = 0;    % Index % initial sample pointer
    numFreqs       = length(freq);                    % length of the frequency bins
    N              = length(rxSig);                   % length of the signal
    nn             = 0:1:N-1;                         % row vector for n

    % Initialization
    sineWave    = exp(-1i*2*pi*freq*nn);
    in_sine     = sineWave(:, 1:N);  % sine "buffer"
    outputMat   = zeros(numFreqs, N);
    delayBuffer = zeros(numFreqs, D);
    sineBuffer  = zeros(numFreqs,D);            % For simplicity (delay buffer for sine)
    outVal      = 0;
    outVec      = [];
    while ( sidx < N)
        sidx = sidx + 1;
        % New samples
        newSamp    = in_sine(:,sidx) * rxSig(sidx) ;
        if D > 0
            buffPtr    = mod(sidx-1, D) + 1;
            buffer_out = delayBuffer(:,buffPtr) .* sineBuffer(:,buffPtr); % previous iter
            sineBuffer(:,buffPtr)  = in_sine(:,sidx);
            delayBuffer(:,buffPtr) = rxSig(sidx) ;
        else
            buffer_out = 0;
        end
        % Out
        outVal     = outVal + newSamp - buffer_out;
        %outVal     =  newSamp ;
        %buffPtrAdj = mod(sidx, D) + 1;
        %outV       = outVal.*conj(sineBuffer(:,buffPtrAdj));
        outV = outVal;
        outVec     = [outVec outV];
    end
    outputMat = outVec;
end

end