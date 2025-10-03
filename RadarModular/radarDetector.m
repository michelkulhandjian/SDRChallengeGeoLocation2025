% Detect Radar signals based on SEC
function [outputMat, det] = radarDetector (data, rP, isFilter)

if ~exist('isFilter','var')
    isFilter = false;
end

det.pMS = nan; det.pMSs=nan;

% MF parameters
N               = rP.Np;
if strcmp(rP.mode,'DETECT')
    outputMat   = data.STFT_sparse;
else
    N           = data.WINDOW;
    outputMat   = data.outputMat;
    %sigAll      = data.dataDet;
    %Nr          = size(sigAll,2);
    Nr          = 1;
end

if isFilter % This part should only select radar samples for ROC should be done in full SEC level
    offsetPW        = -2;
    start_p         = 0;
    lag             = 60 ; %300; % windowSIze of moving average

    SEC_vec         = sum(abs(outputMat));
    SEC_vec_        = movAverageFilter ( SEC_vec, lag);
    SEC_vec_        = SEC_vec_(lag+1:end);
    SEC_vec_        = SEC_vec_ - min(SEC_vec_);

    [cpeaks, TEst, PWEst] = findPeaks (  SEC_vec_');
    Nl              = length( SEC_vec);

    if PWEst < TEst
        indPW     = [-offsetPW - fix(PWEst/2):offsetPW + fix(PWEst/2)];
        LPW       = length(indPW);
        indPulses = (cpeaks(1:end))+indPW';
        indPulses = indPulses(:);
        indPulses(indPulses<1) = [];
        indPulses(indPulses>Nl) = [];
        if ~isempty(indPulses)
            outputMat = outputMat(:, indPulses);
            N     = LPW;
        end
    end
    % save parameters
    det.cpeaks = cpeaks;
    det.TEst = TEst;
    det.PWEst = PWEst;
end

numW  = floor(size(outputMat, 2)/N);

rP.n = 1;
%numW = min(numW, 10);
for n = 1: numW
    % Sparse Output
    rP.tau = 0*(n-1)*N;
    cur_outputMat = outputMat(:,(n-1)*N+1:n*N);
    if strcmp(rP.mode,'DETECT')
        rP.n = n;
        [fig, det.dMF, det.pMSs(1,n), det.estfMF(1,n), det.estgMF(1,n) ] = MF_SEC (cur_outputMat, N, false, data, rP);
    else
        %[fig, det.dMF, det.pMS(1,n), det.estfMF(1,n), det.estgMF(1,n) ] = MF_SEC (cur_outputMat, N, false, data, rP);

        %f0 = 0; tau = 0; Ts = 1/rP.fs;
        %nRange_ = [0:N-1]'; g1 = 1e12;
        %sig= exp( 1i* ( 2*pi*f0*Ts*(nRange_+tau) + pi*g1.*(Ts*(nRange_+tau)).^2 ) );
        for jj = 1 : Nr
            %[~, det.dWVHT, det.pWVHT(1,n), det.estfWVHT(1,n), det.estgWVHT(1,n)] = WVHT (sig, data, rP);
            [det.estfWVHT(1,n), det.estgWVHT(1,(n-1)*Nr + jj)] = CHIRPRATE_SEC (cur_outputMat, N, rP);
            data.indP = (n-1)*N+1:n*N;
            %[~, det.dWVHT, det.pWVHT(1,n), det.estfWVHT(1,n), det.estgWVHT(1,(n-1)*Nr + jj)] = WVHT (sigAll((n-1)*N+1:n*N,jj), data, rP);
            %         [~, det.dHOUGH, det.pHOUGH(1,n), det.estfHOUGH(1,n), det.estgHOUGH(1,n)] = STFT_HOUGH (sig, data, rP);
            %         [~, det.dHOUGH, det.pHOUGH(1,n), det.estfHOUGH(1,n), det.estgHOUGH(1,n)] = STFT_HOUGH (sigAll((n-1)*N+1:n*N,1), data, rP);
            %         [~, det.dMFSTFT, det.pMFSTFT(1,n), det.estfMFSTFT(1,n), det.estgMFSTFT(1,n)] = MF_STFT (sig, data, rP);
            %        [~, det.dMFSTFT, det.pMFSTFT(1,n), det.estfMFSTFT(1,n), det.estgWVHT(1,(n-1)*Nr + jj)] = MF_STFT (sigAll((n-1)*N+1:n*N,1), data, rP);
            %det.estgMF = det.estgWVHT;
            det.estgMF = det.estgWVHT;
        end
    end
end

end
% End of radarDetector

% Moving Average Filter
function [y] = movAverageFilter ( x, windowSize)

PLOT_FIG = false;
%windowSize = 300;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,x);
if PLOT_FIG
    figure; plot(y)
end
end

% chirpRate via SEC output method
function [estf, estg] = CHIRPRATE_SEC (SECOut, N, rP)

PLOT_FIGS = false;
% Initialization of parameters
percent_peak = 0.5;
overlapSmpls = 4;
margin       = 16;

maxSEC = max(max(SECOut));
minSEC = min(min(SECOut));
SEC_m  = zeros(size(SECOut));
[SEC_max]  = SparseP (SECOut', 1, N)';
if sum(find(SEC_max(end-overlapSmpls:end,:))) < sum(find(SEC_max(1:overlapSmpls,:)))
    SEC_max(end-overlapSmpls:end, :) = 0;
else
    SEC_max(1:overlapSmpls,:) = 0;
end
SEC_max = SEC_max(:,margin+1:end-margin);
threshold = minSEC + (maxSEC-minSEC)*percent_peak;
inx  = SEC_max>threshold;
SEC_m(inx)   = SEC_max(inx) ;
%ff = sum(SEC_m);
%[~, time_inx] = find( sum(inx)==1);
[f_inx, time_inx] = find(inx);
%b1 = time_inx\f_inx
%f_inx\time_inx
if ~isempty(time_inx)
    xx = (time_inx-1*time_inx(1))/61.44;
    yy = (f_inx-1*f_inx(1));
    b1 = (xx\yy)*0.7335;

    if b1 < 0.20
        b1 = 1;
    end
else
    b1 = 1;
end
estf = 0;
estg = b1*1e12;

%x_ = mean(xx);
%y_ = mean(yy);
%b2 = sum((xx-x_).*(yy-y_))/sum((xx-x_).^2);

if PLOT_FIGS
    figure; surf(abs(SEC_max))
    figure; plot(time_inx,f_inx )
    figure; surf(abs(SEC_m))
    figure; plot(xx,yy ) ; hold on; plot(xx,b1*xx )
end


end % chirpRate via SEC

% MF via SEC2 output method
function [fig, detectRadar, peakV, estf, estg] = MF_SEC (SECOut, N, sparse_flag, data, rP)

PLOT_FIGS = false;
% Initialization of parameters
D               = data.D;
currf           = data.freq;
totalFreqBins   = data.totalFreqBins;
num_sparse_elem = data.num_sparse_elem;

fig = 0;
tau = rP.tau;
Ts = 1/rP.fs;
gRange = rP.gRange; fRange = rP.fRange;
nRange = [0:N-1];

if isfield(rP,'searchGOnly')
    fRange = rP.fRange1;
end

MF = zeros(length(fRange), length(gRange));
for fInd = 1: length(fRange)
    % grab the f0 initial frequency value for this loop
    f0 = fRange(fInd);
    for gInd = 1: length(gRange)
        % grab the g chip rate value for this loop
        g1 = gRange(gInd);

        if rP.n == 1
            nRange_ = [0:N-1+D];  %( (n-1)*N+1:n*N+start_p);
            sig= exp( 1i* ( 2*pi*f0*Ts*(nRange_+tau) + pi*g1.*(Ts*(nRange_+tau)).^2 ) );
            [curOut_] = SEC_ALG(sig, D, currf, false);
            curOut = curOut_(:,1:end-D);
        else
            dr = 1+mod(rP.n-1, rP.numWperPW);
            nRange_ = [0:N-1]+(dr-1)*N*1;  %( (n-1)*N+1:n*N+start_p);
            sig= exp( 1i* ( 2*pi*f0*Ts*(nRange_+tau) + pi*g1.*(Ts*(nRange_+tau)).^2 ) );
            if dr == 1
                nR_ = [N-D:N-1]+(rP.numWperPW-1)*N*1;
            else
                nR_ = [-D:-1]+(dr-1)*N*1;
            end
            sig_= exp( 1i* ( 2*pi*f0*Ts*(nR_+tau) + pi*g1.*(Ts*(nR_+tau)).^2 ) );
            sig = [sig_ sig];
            [curOut_] = SEC_ALG(sig, D, currf, false);
            curOut = curOut_(:,D+1:end);
        end
        if sparse_flag
            [curOut]  = SparseP (curOut, num_sparse_elem, totalFreqBins);
        end
        SECOut_vec         = sum(abs(SECOut));
        curOut_vec         = sum(abs(curOut));
        conv_SEC = conv(SECOut_vec , fliplr(curOut_vec) );
        MF(fInd, gInd) = max( abs(conv_SEC)) ;
    end
end
if PLOT_FIGS
    if sparse_flag
        titleName = 'MF-SEC-Sparse';
    else
        titleName = 'MF-SEC';
    end
    xlabelName = '$\tilde{g}$ (GHz)';
    ylabelName = '$\tilde{f}$ (MHz)';  zlabelName = 'Normalized MF-SEC Magnitude';
    fig = plotFigs ((gRange)/1e9, (fRange)/1e6, rescale(abs(MF)), titleName, xlabelName, ylabelName, zlabelName );
end

[detectRadar, peakV, estf, estg] = detectDec (MF, gRange, fRange) ;

end % MF via SEC

% Detection decision
function [detectRadar, peakV, estf, estg] = detectDec ( zz, xRange, yRange)

detectRadar = 0; peakV = 0;
peakV = max(max(abs(zz) ));
if size (zz,1) ==  1
    [peakV1,indCol] = max(zz);
    estg = xRange(indCol);
    estf = yRange(1);
elseif size (zz,2) ==  1
    [peakV1,indRow] = max(zz);
    estg = xRange(1);
    estf = yRange(indRow);
else
    [colM,indRow] = max(abs(zz));
    [peakV1,indCol] = max(colM);
    estg = xRange(indCol);
    estf = yRange(indRow(indCol));
end
end % Detection decision

% Plots the figures
function [fig] = plotFigs (xx, yy, zz, titleName, xlabelName, ylabelName, zlabelName )

fig= figure;
%fig = figure('visible','off');
surf(xx,yy,zz,...
    'FaceAlpha',0.9,'EdgeColor','black')
title(titleName,'interpreter','latex','fontsize',18');
xlabel(xlabelName,'interpreter','latex','fontsize',14);
ylabel(ylabelName,'interpreter','latex','fontsize',14);
zlabel(zlabelName,'interpreter','latex','fontsize',14')
set(gca,'TickLabelInterpreter', 'latex');
grid on; set(gcf,'color','w');
%fig = gcf;
end % plotting figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WVHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WVHT Method
function [fig, detectRadar, peakV, estf, estg] = WVHT (s, data, rP)

PLOT_FIGS   = false;
mode        = 1;
fig         = 0;
tau         = rP.tau;
Ts          = 1/rP.fs;      % sample interval
tau         = 0;            % tau
orderSG     = 3;
framelen    = 11;
orderWD     = 4;
[N]         = size(s,1);
if mod(N,2) == 1
    N = N-1;
end
kRange      = -N+1:1:N-1;  % size 2*N - 1
nRange      = [0:N-1]';
NCemp       = 1;
NMF         = 1;
sigAll      = data.dataDet;
[~, Nr]     = size(sigAll);

if (mode == 2)
    NCemp =Nr;
elseif (mode == 3)
    NMF = Nr;
end

gRange = rP.gRange; fRange = rP.fRange;
%nRange = [0:N-1];

if isfield(rP,'searchGOnly')
    fRange = rP.fRange1;
end

MF = zeros(length(fRange), length(gRange));

for mm = 1 : NMF
    % for loop for adding MF matrix
    s = sigAll(data.indP,mm);
    %sR = wdenoise(real(s),orderWD);
    %sC = wdenoise(imag(s),orderWD);
    %s = complex(sR, sC);
    s = sgolayfilt(s,orderSG,framelen);
    s = s(1:N);
    Cemp = funXCorrEmp(s);

    % for loop for adding Cemp
    for jj = 2 : NCemp
        s = sigAll(data.indP,jj);
        s = sgolayfilt(s,orderSG,framelen);
        s = s(1:N);
        % calculate empirical cross-correlation function
        [Cemp]    = Cemp + funXCorrEmp(s)/NCemp;
    end

    for fInd = 1: length(fRange)
        % grab the f0 initial frequency value for this loop
        f0 = fRange(fInd);
        for gInd = 1: length(gRange)
            % grab the g chip rate value for this loop
            g1 = gRange(gInd);

            FF = exp(1i*(f0+g1*(nRange-tau)*Ts)*4*pi*kRange*Ts);
            MF(fInd, gInd) =  MF(fInd, gInd) + WVHTsum3(Cemp,FF, N );
        end
    end
    if PLOT_FIGS
        xlabelName = '$\tilde{g}$ (GHz)';
        ylabelName = '$\tilde{f}$ (MHz)';  zlabelName = 'Normalized MF-SEC Magnitude';
        fig = plotFigs ((gRange)/1e9, (fRange)/1e6, rescale(abs(MF)), titleName, xlabelName, ylabelName, zlabelName );
    end
end

[detectRadar, peakV, estf, estg] = detectDec (MF, gRange, fRange) ;

end % WVHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Emperical Cross-correlation
function [Cemp] = funXCorrEmp( s, nsparse)

if ~exist('nsparse','var')
    % third parameter does not exist, so default it to something
    nsparse = 0;
end

% length of the signal
Nsig = length(s);

kRange = -Nsig+1:1:Nsig-1;
nRange = [0:Nsig-1]';
LenKR = length(kRange);

for kInd = 1:LenKR
    xPlus = circshift(s,-(kInd));
    xMinus = circshift(s,(kInd-LenKR-1));
    Cemp(:,kInd) = xPlus(1:Nsig) .* conj(xMinus(1:Nsig));
end

if nsparse > 0
    Cemp1 = zeros(size(Cemp));
    indR = N:nsparse:2*N-1;
    indL = N-nsparse:-nsparse:1;
    Cemp1(:,indR) = Cemp(:,indR);
    Cemp1(:,indL) = Cemp(:,indL);
    Cemp = Cemp1;
end
% plot real parts to compare
%figure; imagesc(real(Cemp)); title('Real Part - Empirical 1')
end % end of funXCorrEmp

% computing sum of WVHT 3
function [Z] = WVHTsum3 ( Cemp1, Cemp2, N)

%FF = exp(-1i*(ff+gg*(nRange-tau)*Ts)*4*pi*kRange*Ts  );
CF = Cemp1.*conj(Cemp2);

Z = 0;
for n = 0 : N/2-1
    Z = Z + abs(sum ( CF(n+1, N-n:N+n) ) +  sum ( CF(n+1+N/2, N/2+1+n:3/2*N-n-1) ));
end
end % end of WVHTsum3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MF via STFT method
function [fig, detectRadar, peakV, estf, estg] = MF_STFT (s, data, rP)

N  = size(s,1);
% if mod(N,2) == 1
%     N = N-1;
%     s = s(1:N);
% end
orderSG  = 3;
framelen = 11;
orderWD  = 4;
%s = wdenoise(s,orderWD);
sR = wdenoise(real(s),orderWD);
sC = wdenoise(imag(s),orderWD);
%s = complex(sR, sC);
s = sgolayfilt(s,orderSG,framelen);

PLOT_FIGS = false;
% Initialization of parameters
NFFT = 16;
NOVERLAP = 8;
WINDOW = 16;
[~,Ftx,Ttx,PPtx] = spectrogram(s,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
PPtxLog = 10*log10(PPtx);
PPtxLog = PPtxLog/norm(PPtxLog, 'fro');

fig = 0;
tau = rP.tau;
Ts = 1/rP.fs;
gRange = rP.gRange; fRange = rP.fRange;
nRange = [0:N-1];

if isfield(rP,'searchGOnly')
    fRange = rP.fRange1;
end

MF = zeros(length(fRange), length(gRange));
for fInd = 1: length(fRange)
    % grab the f0 initial frequency value for this loop
    f0 = fRange(fInd);
    for gInd = 1: length(gRange)
        % grab the g chip rate value for this loop
        g1 = gRange(gInd);

        sig= exp( 1i* ( 2*pi*f0*Ts*(nRange+tau) + pi*g1.*(Ts*(nRange+tau)).^2 ) );
        [~,curFtx,curTtx,curPPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
        curPPtxLog = 10*log10(curPPtx);
        curPPtxLog = curPPtxLog/norm(curPPtxLog, 'fro');
        [norm(s-sig.'), norm(PPtxLog - curPPtxLog )];
        %MF(fInd, gInd) = MF2D (PPtxLog, curPPtxLog, 6, 6);
        MF(fInd, gInd) = sig*conj(s);
    end
end
if PLOT_FIGS
    xlabelName = '$\tilde{g}$ (GHz)';
    ylabelName = '$\tilde{f}$ (MHz)';  zlabelName = 'Normalized MF-SEC Magnitude';
    fig = plotFigs ((gRange)/1e9, (fRange)/1e6, rescale(abs(MF)), titleName, xlabelName, ylabelName, zlabelName );
end

[detectRadar, peakV, estf, estg] = detectDec (MF, gRange, fRange) ;

end % MF via STFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peakV] = MF2D (CC, ZZ, xShift, yShift)

[Ny, Nx]  = size (CC);
ZZ = conj(ZZ);
MF   = zeros(yShift*2+1, xShift*2+1);
tempC  = zeros (size(CC));
tempZ  = zeros(size(ZZ));
yRange = -yShift : yShift;
xRange = -xShift : xShift;
for indY = 1 : length(yRange)
    yy = yRange(indY);
    for indX = 1 : length(xRange)
        xx = xRange(indX);
        % C = tempC;
        Z = tempZ;
        %% y
        indYZ1  = (1)*(sign(1-yy-0.5)+1)/2+(1+yy)*(sign(yy-0.5)+1)/2;
        indYZ2  = (Ny+yy)*(sign(1-yy-0.5)+1)/2+(Ny)*(sign(yy-0.5)+1)/2;
        indYZZ1 = (1-yy)*(sign(1-yy-0.5)+1)/2+(1)*(sign(yy-0.5)+1)/2;
        indYZZ2 = (Ny)*(sign(1-yy-0.5)+1)/2+(Ny-yy)*(sign(yy-0.5)+1)/2;
        %% x
        indXZ1  = (1)*(sign(1-xx-0.5)+1)/2+(1+xx)*(sign(xx-0.5)+1)/2;
        indXZ2  = (Nx+xx)*(sign(1-xx-0.5)+1)/2+(Nx)*(sign(xx-0.5)+1)/2;
        indXZZ1 = (1-xx)*(sign(1-xx-0.5)+1)/2+(1)*(sign(xx-0.5)+1)/2;
        indXZZ2 = (Nx)*(sign(1-xx-0.5)+1)/2+(Nx-xx)*(sign(xx-0.5)+1)/2;
        Z (indYZ1:indYZ2 , indXZ1:indXZ2 ) =  ZZ (indYZZ1:indYZZ2, indXZZ1:indXZZ2);
        MF(indY, indX) = sum(sum(CC.*Z));
    end
end

peakV = max(max(MF));

end
% MF 2-D



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STFT-HOUGH method
function [fig, detectRadar, peakV, estf, estg] = STFT_HOUGH (s, data, rP)

N           = size(s,1);
if mod(N,2) == 1
    N = N-1;
    s = s(1:N);
end

PLOT_FIGS = false;
% Initialization of parameters
NFFT = 16;
NOVERLAP = 8;
WINDOW = 16;
[~,Ftx,Ttx,PPtx] = spectrogram(s,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
PPtxLog = 10*log10(PPtx);
[P, lines] = lineDetection ( PPtxLog);

fig = 0;
tau = rP.tau;
Ts = 1/rP.fs;
gRange = rP.gRange; fRange = rP.fRange;
nRange = [0:N-1];

if isfield(rP,'searchGOnly')
    fRange = rP.fRange1;
end

MF = zeros(length(fRange), length(gRange));
for fInd = 1: length(fRange)
    % grab the f0 initial frequency value for this loop
    f0 = fRange(fInd);
    for gInd = 1: length(gRange)
        % grab the g chip rate value for this loop
        g1 = gRange(gInd);

        sig= exp( 1i* ( 2*pi*f0*Ts*(nRange+tau) + pi*g1.*(Ts*(nRange+tau)).^2 ) );
        [~,curFtx,curTtx,curPPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
        curPPtxLog = 10*log10(curPPtx);
        [curP, curlines] = lineDetection ( curPPtxLog);
        MF(fInd, gInd) = ErrorLineEst (lines, curlines) ;
    end
end
if PLOT_FIGS
    xlabelName = '$\tilde{g}$ (GHz)';
    ylabelName = '$\tilde{f}$ (MHz)';  zlabelName = 'Normalized MF-SEC Magnitude';
    fig = plotFigs ((gRange)/1e9, (fRange)/1e6, rescale(abs(MF)), titleName, xlabelName, ylabelName, zlabelName );
end

[detectRadar, peakV, estf, estg] = detectDec (MF, gRange, fRange) ;

end % STFT-HOUGH


% HT -Line Detection
function [P, lines] = lineDetection ( ZZ)

example_flag = false;
plot_figs    = false;

if example_flag
    % create a white image of size 200x200
    im = uint8 (zeros(201, 201)) + 255;
    d = [0,0; 200,0;0,200;200,200; 100,100;100,30];
    x1 = 1;
    RGB = im;
    while (x1 < 6)
        n = 6-x1;
        for i = 1 : n
            RGB = insertShape (RGB, 'Line', [d(x1,1)+1, d(x1,2)+1, d(x1+i,1)+1, d(x1+i,2)+1],...
                'Color', [0 0 0], 'Opacity', 1);
        end
        x1 = x1 + 1;
    end
    imshow(imbinarize(rgb2gray(RGB)));

    % Converting Image to BW
    Img = RGB;
    Img2 = rgb2gray (Img);
    Img2 = imbinarize(Img2);
    BW = imcomplement(Img2);
    imshow(BW);
    ZZ = BW;
end

[H,theta,rho] = hough(ZZ);
if plot_figs
    figure; imshow(imadjust(rescale(H)),'XData',theta,'YData',rho,...
        'InitialMagnification','fit');
    title('Hough transform of STFT_sparse');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    colormap(gca,hot);
end

% Applying threshold to find points on lines
threshold = 0.5;
P = houghpeaks(H, 15, 'threshold', ceil(threshold*max(H(:))));
x = theta(P(:,2));
y = rho(P(:,1));
%lines = houghlines (ZZ, theta, rho, P, 'FillGap', 20, 'MinLength', 7);
lines = houghlines(ZZ,theta,rho,P) ;
if plot_figs
    figure, hold on;
    for k = 1: length(lines)
        xy = [lines(k).point1 ; lines(k).point2];
        plot(xy(:,1), xy(:,2), 'LineWidth', 2 );
        % Plot beginnings and eds of lines
        plot( xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'green');
        plot( xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'red');

    end
end
end % HT -Line Detection

% HT -Line Detection
function [err] = ErrorLineEst ( lines, curlines)

err = 0;
if length(lines) >0
    if length(curlines) == 0
        err = 1e5;
    else
        for k = 1: min(length(lines), length(curlines) )
            err = err+ norm(lines(k).point1-curlines(k).point1) +  norm(lines(k).point2-curlines(k).point2);
        end
    end
end
end

