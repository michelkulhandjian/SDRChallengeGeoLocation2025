
function[IQout,  bufLen, fc, packet_len] = gen_signal ( mode, bw, fs, Ts, rP) 

%params
txfc  = rP.fc/1e9;        % Tx Center frequency (GHz) e.g., 2.5 GHz  
period = Ts;              %  Period (s) e.g.,  0.1024
ModL = mode ;             % Modulation
fs  = fs/1e6;             % Sampling rate (Hz) e.g., 40 MHz
bw =  bw/1e6;             % Bandwidth (Hz) e.g.,  10 MHz
fc = rP.fcLoc ;           % Local Center frequency (MHz) e.g., -15 MHz
pwr = rP.pwr;             % Power  Peak power (dB) e.g., -40 dB
beta = rP.beta;           % Excess Bandwidth e.g.,   0.05
NFFT =  rP.NFFT;          % Size of FFT  e.g., -32768
out_cmplx = rP.out_cmplx; % Output IQ e.g., CMPLX' or 'REAL'
PLOT_FIG  = true;       


NTswitch = 20;
Tswitch = zeros(1, NTswitch); % April 14, 2025, Kazi made it Dynamic length based on Tswitch_cell
Tswitch_len = 0;
%signals = 0;

% Set bps - bits / symbol and the constelation table

if strcmp(ModL, 'BPSK')
    mod_type = 'QAM';
    bps = 1;
    ctable = qammod((0:1),2,"UnitAveragePower",true);    
elseif strcmp(ModL, 'QPSK')
    mod_type = 'QAM';
    bps = 2;
    ctable = qammod((0:3),4,"UnitAveragePower",true);   
elseif strcmp(ModL, '8PSK')
    mod_type = 'QAM';
    bps = 3;
    ctable = pskmod((0:7),8,0,'gray',PlotConstellation=true,InputType='integer');
elseif strcmp(ModL, 'D8PSK')
    mod_type = 'QAM';
    bps = 3;
    ctable = dpskmod((0:7),8);
elseif strcmp(ModL, '16PSK')
    mod_type = 'QAM';
    bps = 4;
    ctable = pskmod((0:15),16,0,'gray',PlotConstellation=true,InputType='integer');
elseif strcmp(ModL, 'QAM16')
    mod_type = 'QAM';
    bps = 4;
    ctable = qammod((0:15),16,"UnitAveragePower",true);
 elseif strcmp(ModL, 'QAM32')   
    mod_type = 'QAM';
    bps = 5;    
    ctable = qammod((0:31),32,"UnitAveragePower",true);
 elseif strcmp(ModL, 'QAM64')
    mod_type = 'QAM';
    bps = 6;
    ctable = qammod((0:63),64,"UnitAveragePower",true);
elseif strcmp(ModL, 'QAM256')
    mod_type = 'QAM';
    bps = 8;
    ctable = qammod((0:255),256,"UnitAveragePower");
elseif strcmp(ModL, 'NXDN48')
    mod_type = 'FSK';
    bps = 2;
    bd = 2400;      % Baud rate
    fmod = 350*[1, 3, -1,-3];
elseif strcmp(ModL, 'NXDN96')
    mod_type = 'FSK';
    bps = 2;
    bd = 4800;      % Baud-rate
    fmod = 800*[1, 3, -1,-3];
elseif strcmp(ModL, 'TONE')
    mod_type = 'TONE';
    bps = 1;
elseif strcmp(ModL, 'PAM2')
    mod_type = 'PAM';
    bps = 1;
    ctable = pammod((0:1),2);
elseif strcmp(ModL, 'PAM4')
    mod_type = 'PAM';
    bps = 2;
    ctable = pammod((0:3),4);
elseif strcmp(ModL, 'PAM8')
    mod_type = 'PAM';
    bps = 3;
    ctable = pammod((0:7),8);
elseif strcmp(ModL, 'PAM16')
    mod_type = 'PAM';
    bps = 4;
    ctable = pammod((0:15),16);
end

if strcmp(mod_type, 'QAM') || strcmp(mod_type, 'PAM')
    if ismissing(fs)
        %os is the over sampling parameter
        if txfc/bw*1e3 > 10
            os = 4;
        else
            os = 8;
        end
        fs = bw * os;
    else
        os = fs / bw;
    end
    packet_len = fs * (period * 1e6);
    symlen = floor(packet_len / os);
    bufLen = round(symlen * os);

elseif strcmp(mod_type, 'FSK')
    if ismissing(fs)
        os = 4;
        fs = bd*os*1e-6;
    else
        os = fs/bd*1e6;
    end
    % G = gcd(fres*1e6, fs*1e6);
    % P = fres*1e6 / G;
    % Q = fs*1e6 / G;
    symlen = floor(period * bd);
    packet_len = fs * (period * 1e6);
    % packet_len_res = packet_len * fres / fs;
elseif strcmp(mod_type, 'TONE')
    os = 1; %os is the over sampling parameter
    % P=1;
    % Q=1;
    packet_len = fs * (period * 1e6);
    % packet_len_res = packet_len;
    bufLen = packet_len;
    symlen = bufLen;
end

% Set up the first pulse shaping filter
if ~(strcmp(mod_type, 'FSK'))
    rspan = 60;
    hr = firls(round(os*rspan)-1, [0, (1-beta/2)*bw/fs,(1+beta/2)*bw/fs, 1], [1,1,0,0]);
else
    rspan = 8;
    hr = firls(round(os*rspan)-1, [0, (1-beta/2)*bw/fs,(1+beta/2)*bw/fs, 1], [1,1,0,0]);
    
end

% Set up the second pulse shaping filter
tspan = 10;
if ~(strcmp(mod_type, 'FSK'))
    ht = rcosdesign(beta,2*tspan,round(os), "normal");       % The RRC
else
    ht = rcosdesign(0.8,2*tspan,round(os/2), "normal");       % The RRC    
end    

if (strcmp(mod_type, 'QAM')) || (strcmp(mod_type, 'FSK')) || strcmp(mod_type, 'PAM')
    h = conv(ht, hr);
    h = h / sqrt(h*h');       % Normalize h
    span = rspan + 2*tspan;
else
    h = 1;
    span = 1;
end

% Convert to analitic filter
ha = hilbert(h);

if strcmp(mod_type, 'FSK')
    bufLen = round(symlen * os);
end

if strcmp(out_cmplx, 'Real')
    S_buf = zeros(1, bufLen);   % Input buffers to h
    I_buf = zeros(1, bufLen);
else
    Q_buf = zeros(1, bufLen);   % Input buffers to h
    I_buf = zeros(1, bufLen);
end
IQ_merged = zeros(1, 2*bufLen);
IQ_merged_test = zeros(1, 2*bufLen);
symPtr = 0;
wr_ptr = 0;
IQout = zeros(1, bufLen);
dbuf = zeros(symlen, bps);
init_phase = 0;     % For FSK

% Fill the data buffers with data
% Set up the generation of a string of random data, uniformly distributed 0 and 1
rng('default');
if strcmp(mod_type, 'QAM' ) || strcmp(mod_type, 'PAM')

    %Generate all the symbols :symlen
    
    din = randi([0,1], [symlen, bps]);
    dbuf(1:symlen,:) = din;
    % Convert din to symbol
    symin = bin2dec(num2str(din));

    % Get the Real and Imag from the constelation table
    re = real(ctable(symin+1));
    im = imag(ctable(symin+1));

    I_buf(1 : symlen) = re;
    Q_buf(1 : symlen) = im;
    
    if strcmp(out_cmplx, 'Real')
        % Run it through H filter and get os outputs
        s_state = S_buf(mod((-span : -1), bufLen)+1);
        sinqam = re + 1i*im;
        shin = kron([s_state, sinqam], [1, zeros(1, os-1)]);
        filtout = filter(ha, 1, shin);
        % ihout = real(filtout);
        sout = filtout(end-symlen*os+1 : end);
        % Save the outputs in the I and Q buffers
        IQ(1 : symlen*os) = real(sout);
        
    else
        % Run it through H and get os outputs
        q_state = Q_buf(mod((-span : -1), bufLen)+1);
        i_state = I_buf(mod((-span : -1), bufLen)+1);
        qhin = kron([q_state, im], [1, zeros(1, os-1)]);
        ihin = kron([i_state, re], [1, zeros(1, os-1)]);
        qhout = filter(h,1,qhin);
        ihout = filter(h,1,ihin);
        qout = qhout(end-symlen*os+1 : end);
        iout = ihout(end-symlen*os+1 : end);
        
        % Save the outputs in the I and Q buffers
        IQ_merged(1 : 2 : 2*symlen*os) = iout;
        IQ_merged(2 : 2 : 2*symlen*os) = qout;
    end

    if ~strcmp(out_cmplx, 'Real')
        IQ = IQ_merged(1:2:end) + 1i*IQ_merged(2:2:end);
    end


elseif strcmp(mod_type, 'FSK')

    wr_ptr = 0; %writer pointer for symbols

    for symPtr = 1 : symlen*2-2 % Corrected Kazi March 12, 2025  % for symPtr = 1 : symlen 

        %Generate a symbol using uniform distribution 
        din = randi([0,1], [1, bps]);
        dbuf(symPtr,:) = din;
        % Convert din to symbol
        symin = bin2dec(num2str(din));
        t = (0:os-1);
        re = cos(2*pi*fmod(symin+1)*t/(fs*1e6*os) + init_phase);
        init_phase = 2*pi*fmod(symin+1)*os/(fs*1e6*os) + init_phase;
        I_buf(wr_ptr + 1 : wr_ptr + os) = re;        
        
        if ~strcmp(out_cmplx, 'Real')
            % Run it through H filter 
            i1_state = I_buf(round(mod((wr_ptr-length(ha) : wr_ptr-1), packet_len)+1));% March 12, 2025 Mushtaq %i1_state = round(I_buf(mod((wr_ptr-length(ha) : wr_ptr-1), packet_len)+1));
            ihin = [i1_state, re];
            ihout = filter(ha,1,ihin);
            iout = ihout(end-os+1 : end);
            IQ_merged(wr_ptr+1 : 2 : wr_ptr+2*length(iout)) = real(iout);
            IQ_merged(wr_ptr+2 : 2 : wr_ptr+2*length(iout)) = imag(iout);
        else
            % Run it through H filter 
            i1_state = I_buf(mod((wr_ptr-length(h) : wr_ptr-1), packet_len)+1);
            ihin = [i1_state, re];
            ihout = filter(h,1,ihin);
            iout = ihout(end-os+1 : end);
            IQ_merged(wr_ptr+1 : 2 : wr_ptr+2*length(iout)) = iout;
            IQ_merged(wr_ptr+2 : 2 : wr_ptr+2*length(iout)) = 0;
        end

        % Save the outputs in the I and Q buffers

        % Update the pointers
        wr_ptr = wr_ptr + os;

    end
    IQ = IQ_merged(1:2:end) + 1i*IQ_merged(2:2:end);
    
elseif  (strcmp(mod_type, 'TONE'))
    t = (0:(1/fs):period*10e6);
    IQ_before = sin(2*pi*fc*t);
    rspan = 100;
    beginning = (fc - bw/2)/(fs/2);
    myend = (fc + bw/2)/(fs/2);
    hr = firls(round(rspan)-1, [0, .99*beginning,1.01*beginning, .99*myend,1.01*myend, 1], [0,0,1,1,0,0]);
    if ~strcmp(out_cmplx, 'Real')
        % Run it through Hilbert filter 
        % Convert to analitic filter
        ha = hilbert(hr);
        IQ = conv(IQ_before,ha);
    else
        % Run it through H filter 
        IQ = conv(IQ_before,hr);
    end
    IQout = zeros(1, length(IQ));
end
% Make each signals power according to the input power vector
%IQ_pwr = (IQ * IQ'); % Mushtaq 03282025
IQ_pwr = (IQ * IQ')/length(IQ);
% if ~strcmp(ModL, 'TONE')
%     IQpsd = 10*log10(IQ_pwr / (bw*NFFT/fs));
% else
% IQpsd = 10*log10(IQ_pwr/bw*fs);
IQpsd = 10*log10(IQ_pwr);
% end

IQscaled = IQ * 10^((pwr - IQpsd)/20);

IQout(1, :) = IQscaled;


if PLOT_FIG
    figure;

    % subplot(2,1,1);
    plot(real(IQout),'g');
    axis auto
    grid on;
    title('Time Domain Signal Before Tswitch - real');
end

% subplot(2,1,2);
% plot(imag(IQout),'r');
% grid on;
% title('QAM32time domain - imag');

% Use the Tswitch vector in excell input to turn the signal on and off in ms

if Tswitch_len > 0
    Ton = Tswitch(1:2:Tswitch_len);
    Toff = Tswitch(2:2:Tswitch_len);
    %Tgap = [1*fs*1e3:Ton(1)*fs*1e3]; %Tgap = [1:Ton(1)]*fs*1e3;
    Tgap = [1:Ton(1)*fs*1e3]; % Kazi Corrected 02142025
    IQout(:, Tgap) = 0;
    
    % Turn the signal off at the off times
    for it = 1 : Tswitch_len/2 - 1
        Tgap = [Toff(it)*fs*1e3:Ton(it+1)*fs*1e3]; 
        IQout(:,round(Tgap)) = 0;
    end
     Tgap = [Toff(Tswitch_len/2)*fs*1e3 : packet_len]; 
     IQout(:,Tgap) = 0; 
end



if PLOT_FIG
    figure;

    % subplot(2,1,1);
    plot(real(IQout),'g');
    axis auto
    grid on;
    title('Time Domain Signal After Tswitch - real');

    % subplot(2,1,2);
    % plot(imag(IQout),'r');
    % grid on;
    % title('QAM32time domain - imag');
    %Test_Signal_Generation=0;
end