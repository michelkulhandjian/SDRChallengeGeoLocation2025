function [rP] = localize2 ( data, rP)

addpath '..\';

% Parameters
c           = physconst('LightSpeed');
SNR         = [0:4:35]+15; % in dB
PLOT_FIG    = true;

%% Create BS Platform Objects
pPos       = rP.pPos;  
bPos       = rP.bPos;
%[rxSig]    = data.sig;
[L, N]     = size(bPos);

% Visualize the transmitter and receiving nodes
viewTxRxposition (bPos, pPos);

%% Testing
% Verify
if 0
    dist = tDelay*c;
    delta = dist(2:end)'-dist(1); delta = [0;delta]

    n = 4
    (pPos(1,1)-bPos(1,1))*(bPos(1,n)-bPos(1,1)) + (pPos(2,1)-bPos(2,1))*(bPos(2,n)-bPos(2,1)) + (pPos(3,1)-bPos(3,1))*(bPos(3,n)-bPos(3,1)) + delta(n,1)*dist(1)
    (1/2)*((bPos(1,n)-bPos(1,1)).^2 + (bPos(2,n)-bPos(2,1)).^2 + (bPos(3,n)-bPos(3,1)).^2 - delta(n,1).^2)

    n= 2:4;
    A0 = [(bPos(1,n)-bPos(1,1))',  (bPos(2,n)-bPos(2,1))',  (bPos(3,n)-bPos(3,1))',  delta(n,1)];
    th0 = [(pPos(1,1)-bPos(1,1));  (pPos(2,1)-bPos(2,1));  (pPos(3,1)-bPos(3,1));  dist(1)];
    A0*th0
    (1/2)*((bPos(1,n)-bPos(1,1))'.^2 + (bPos(2,n)-bPos(2,1))'.^2 + (bPos(3,n)-bPos(3,1))'.^2 - delta(n,1).^2)


    delta(2,1) + dist(1)
    dist(2)
end


% Testing
% true values
if 0
    A0= [bPos(1,2:end)'-bPos(1,1), bPos(2,2:end)'-bPos(2,1), bPos(3,2:end)'-bPos(3,1), range1'];
    th0 = [pPos(1,1)-bPos(1,1); pPos(2,1)-bPos(2,1); pPos(3,1)-bPos(3,1); tDelay(1)*c];
    b0 = 0.5*[(bPos(1,2:end)'-bPos(1,1)).^2+(bPos(2,2:end)'-bPos(2,1)).^2+(bPos(3,2:end)'-bPos(3,1)).^2-(range1.^2)'];
    [A0*th0 b0 A0*th0-b0]
    tt = pinv(A0)*b0;

    % estimated
    A01= [bPos(1,2:end)'-bPos(1,1), bPos(2,2:end)'-bPos(2,1), bPos(3,2:end)'-bPos(3,1), (Y*c)'];
    b01 = 0.5*[(bPos(1,2:end)'-bPos(1,1)).^2+(bPos(2,2:end)'-bPos(2,1)).^2+(bPos(3,2:end)'-bPos(3,1)).^2-((Y*c).^2)'];
    b011 = 0.5*[-norm(bPos(:,1))^2+sum(bPos(:,2:end).^2,1)'-((Y*c).^2)'];
    tt1 = pinv(A01)*b01;
    tt11 = pinv(A01)*b011;
end


switch rP.mode
    %% Geolocate a position
    case "GeolocatePosition"
        %% Alg1
        if rP.isALg1
            pos00 = tdoaposestf(rP.Y,rP.var,bPos);
            display("RMSE for 'Alg 1");
            rmse(pos00,pPos)
            [pos00]
            rP.pos00 = pos00;
            %fprintf(rP.fileID,['Alg 1 ', int2str(pos00'), ' \n']);
             if ~rP.isEuclidean
                [geopos00] = enuTogeo(pos00, rP);
                [geopos00]
                rP.geopos00 = geopos00;
            end
        end
        tdoa = rP.tdoaData;
        [pos0] = geolocationA (bPos, tdoa', N, c,  0); % Linear Solution LS 1
        [pos1] = geolocationB (bPos, tdoa', N, c,  0); % CTLS
        [pos2] = geolocationB (bPos, tdoa', N, c,  1); %  OSNM-CTLS Solution
        [pos3] = geolocationB (bPos, tdoa', N, c,  2, 0);  % Taylor Series
        [pos4] = geolocationB (bPos, tdoa', N, c,  3);  % Chan's Ho
        [pos5] = geolocationB (bPos, tdoa', N, c,  4);  % Improved Chan's Ho
        %[pos6] = geolocationB (bPos, Y', N, c,  5);  % 2-Step Weighting LS solution
        [pos7] = geolocationA (bPos, tdoa', N, c,  2); % Closed form LS solution
        if rP.isNonLinear
            [pos8] = geolocationA (bPos, [tdoa rP.tdoa32]', N, c,  6); % Nonlinear solution
            display("RMSE for 'Nonlinear");
            rmse(pos8,pPos)
            [pos8]
            rP.pos8 = pos8;

            if ~rP.isEuclidean
                [geopos8] = enuTogeo(pos8, rP);
                [geopos8]
                rP.geopos8 = geopos8;
            end
        end

        display("RMSE for 'Traditional LS', 'CTLS',  'CTLS Iterative', 'Taylor Series', 'CHans Ho', 'Imp. Chans Ho', 'Closed form LS'")
        [ rmse(pos0,pPos), rmse(pos1,pPos), rmse(pos2,pPos), rmse(pos3,pPos), rmse(pos4,pPos), rmse(pos5,pPos), rmse(pos7,pPos)]
        [pos0, pos1, pos2, pos3, pos4,  pos7]
        rP.pos0 = pos0; rP.pos1= pos1; rP.pos2 = pos2; rP.pos3 = pos3; rP.pos4 = pos4;
        rP.pos5 = pos5; rP.pos7 = pos7;
        if ~rP.isEuclidean
            [geopos0] = enuTogeo(pos0, rP); [geopos1] = enuTogeo(pos1, rP); [geopos2] = enuTogeo(pos2, rP);
            [geopos3] = enuTogeo(pos3, rP); [geopos4] = enuTogeo(pos4, rP); [geopos5] = enuTogeo(pos5, rP);
            [geopos7] = enuTogeo(pos7, rP);
            [geopos0, geopos1, geopos2, geopos3, geopos4, geopos5, geopos7]
            rP.geopos0 = geopos0; rP.geopos1= geopos1; rP.geopos2 = geopos2; rP.geopos3 = geopos3; rP.geopos4 = geopos4;
            rP.geopos5 = geopos5; rP.geopos7 = geopos7;
        end

        %% Geolocate multiple positions
    case "GeolocateMultiplePosition"
        adfd

        %% Simulation Geolocation Algorithms.
    case "GeolocateSimulation"

        NN = 1000;  % number of iterations

        RMSE_A = zeros(size(SNR));
        RMSE_B = zeros(size(SNR));
        RMSE_C = zeros(size(SNR));
        RMSE_D = zeros(size(SNR));
        RMSE_E = zeros(size(SNR));
        RMSE_F = zeros(size(SNR));
        RMSE_G = zeros(size(SNR));
        RMSE_I = zeros(size(SNR));
        RMSE_J = zeros(size(SNR));
        RMSE_K = zeros(size(SNR));

        POS00 = zeros(length(pPos),NN,length(SNR) );
        POS0 = zeros(length(pPos),NN,length(SNR) );
        POS1 = zeros(length(pPos),NN,length(SNR) );
        POS2 = zeros(length(pPos),NN,length(SNR) );
        POS3 = zeros(length(pPos),NN,length(SNR) );
        POS4 = zeros(length(pPos),NN,length(SNR) );
        POS5 = zeros(length(pPos),NN,length(SNR) );
        POS6 = zeros(length(pPos),NN,length(SNR) );
        %POS7 = zeros(length(pos0),NN,length(SNR) );
        POS8 = zeros(length(pPos),NN,length(SNR) );

        tdoa = rP.tdoaData; 
        for ind = 1: length(SNR)
            for nn = 1: NN

                Ps = norm(tdoa);
                Pn = 1*Ps*10^(-SNR(ind)/20);
                tdoa_ = tdoa +  1*Pn*randn(1,N-1);  % adding AWGN to TDOAs as errors

                [pos0] = geolocationA (bPos, tdoa_', N, c,  2);
                if rP.isALg1
                    pos00 = tdoaposestf(tdoa_,rP.var,bPos);
                    % RMSE of the TDOA position estimate
                    RMSE00 = rmse(pos00,pPos);
                    %disp(['RMS Localization error = ', num2str(RMSE), ' meters.'])
                    POS00(:,nn,ind) = pos00;
                    RMSE_A(ind) = RMSE_A(ind)+ RMSE00/NN;
                end
                [pos1] = geolocationB (bPos, tdoa_', N, c,  0);
                [pos2] = geolocationB (bPos, tdoa_', N, c,  1);
                [pos3] = geolocationB (bPos, tdoa_', N, c,  2, 0);
                [pos4] = geolocationB (bPos, tdoa_', N, c,  3);
                [pos5] = geolocationB (bPos, tdoa_', N, c,  4);  % Improved Chan's Ho
                % [pos6] = geolocationB (bPos, Y_', N, c,  5);  % 2-Step Weighting LS solution
                [pos7] = geolocationA (bPos, tdoa_', N, c,  2); % Closed form LS solution


                RMSE0 = rmse(pos0,pPos);
                RMSE1 = rmse(pos1,pPos);
                RMSE2 = rmse(pos2,pPos);
                RMSE3 = rmse(pos3,pPos);
                RMSE4 = rmse(pos4,pPos);
                RMSE5 = rmse(pos5,pPos);
                %RMSE6 = rmse(pos6,pPos);
                RMSE7 = rmse(pos7,pPos);

                if rP.isNonLinear
                    [pos8] = geolocationA (bPos, [tdoa_ rP.tdoa32+1*Pn*randn]', N, c,  6); % Nonlinear solution
                    RMSE8 = rmse(pos8,pPos);
                    POS8(:,nn,ind) = pos8;
                    RMSE_K(ind) = RMSE_K(ind) + RMSE8/NN;
                end

                POS0(:,nn,ind) = pos0;
                POS1(:,nn,ind) = pos1;
                POS2(:,nn,ind) = pos2;
                POS3(:,nn,ind) = pos3;
                POS4(:,nn,ind) = pos4;
                POS5(:,nn,ind) = pos5;
                %POS6(:,nn,ind) = pos6;
                %OS7(:,nn,ind) = pos7;
                
                RMSE_B(ind) = RMSE_B(ind) + RMSE0/NN;
                RMSE_C(ind) = RMSE_C(ind) + RMSE1/NN; RMSE_D(ind) = RMSE_D(ind) + RMSE2/NN;
                RMSE_E(ind) = RMSE_E(ind) + RMSE3/NN; RMSE_F(ind) = RMSE_F(ind) + RMSE4/NN;
                RMSE_G(ind) = RMSE_G(ind) + RMSE5/NN;
                RMSE_J(ind) = RMSE_J(ind) + RMSE7/NN;
            end
        end

        indSNR = 5;

        figure;
        if rP.isALg1
            scatter(POS00(1,:,indSNR), POS00(2,:,indSNR), 'LineWidth',1);
        end

        hold on;

        scatter(POS0(1,:,indSNR), POS0(2,:,indSNR),  'LineWidth',1);
        scatter(POS1(1,:,indSNR), POS1(2,:,indSNR),  'LineWidth',1);
        scatter(POS2(1,:,indSNR), POS2(2,:,indSNR),  'LineWidth',1);
        scatter(POS3(1,:,indSNR), POS3(2,:,indSNR),  'LineWidth',1);
        scatter(POS4(1,:,indSNR), POS4(2,:,indSNR),  'LineWidth',1);
        scatter(POS5(1,:,indSNR), POS5(2,:,indSNR),  'LineWidth',1);
        scatter(POS6(1,:,indSNR), POS6(2,:,indSNR),  'LineWidth',1);
        %scatter(POS7(1,:,indSNR), POS7(2,:,indSNR),  'LineWidth',1);
        if rP.isNonLinear
            scatter(POS8(1,:,indSNR), POS8(2,:,indSNR),  'LineWidth',1);
        end
        plot(pPos(1), pPos(2), 'r*','LineWidth',1,'MarkerSize',10); grid on
        
        if rP.isALg1
            legend(  'Alg 1', 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS', 'Target');
             if rP.isNonLinear
                 legend(  'Alg 1', 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS', 'Nonlinear', 'Target');
             end
        else
            legend( 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS', 'Target');
            if rP.isNonLinear
                legend( 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS', 'Nonlinear', 'Target');
            end
        end


        CRLB  = zeros(3,1);
        Ytotal = 0;
        Xtotal = 0;
        XYtotal = 0;

        for i = 1 : N
            ysq = (pPos(2)-bPos(2,i))^2;
            xsq = (pPos(1)-bPos(1,i))^2;
            xysq = (pPos(1)-bPos(1,i))*(pPos(2)-bPos(2,i));
            totsq = ysq + xsq;
            Ytotal = Ytotal + ysq/totsq;
            Xtotal = Xtotal + xsq/totsq;
            XYtotal = XYtotal + xysq/totsq;

        end
        CRLB(1) = Ytotal/(Xtotal*Ytotal - XYtotal^2);
        CRLB(2) = Xtotal/(Xtotal*Ytotal - XYtotal^2);

        R = 0.5;
        [xc, yc] = drawCircle(pPos, R);
        plot(xc, yc, 'k','LineWidth',3,'MarkerSize',10);
        atitle = sprintf( 'CRLB(x) = %.2f, CRLB(y) = %.2f, SNR =  %.2f dB \n',CRLB(1), CRLB(2), SNR(indSNR) );
        title(atitle)


        %figure; plot(SNR, RMSE_A); hold on; plot(SNR, RMSE_B); plot(SNR, RMSE_C); plot(SNR, RMSE_D);  grid on;
        figure;
        if rP.isALg1
            plot(SNR, RMSE_A, '-*', 'LineWidth',1.5);
        end
        hold on; plot(SNR, RMSE_B, '--', 'LineWidth',1.5); plot(SNR, RMSE_C, '-v', 'LineWidth',1.5); plot(SNR, RMSE_D, '-x', 'LineWidth',1.5);
        plot(SNR, RMSE_E, '-o', 'LineWidth',1.5); plot(SNR, RMSE_F, '-', 'LineWidth',1.5);  plot(SNR, RMSE_G, '-s', 'LineWidth',1.5);  plot(SNR, RMSE_J, '-x', 'LineWidth',1.5);  grid on;
        if rP.isNonLinear
            plot(SNR, RMSE_K, '-v', 'LineWidth',1.5); % plot(SNR, RMSE_I, 'k', 'LineWidth',1.5);
        end
        xlabel('SNR (dB)') ; ylabel('RMSE (meters)'); 
       if rP.isALg1
            legend(  'Alg 1', 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS');
            if rP.isNonLinear
                legend(  'Alg 1', 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS', 'Nonlinear');
            end
       else
            legend( 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg',  'Closed form LS');
            if rP.isNonLinear
                legend( 'Traditional LS',  'CTLS',  'CTLS Iterative', 'Taylor Series', 'Chans Ho' , 'Imp. Chans Alg', 'Closed form LS', 'Nonlinear');
            end
       end
        %atitle = sprintf( 'Localisation Algorithms No. Rx: %d, GDOP = %.2f  \n', numRadar, GDOP);
        atitle = sprintf( 'Localisation Algorithms No. Rx: %d \n', N);
        title(atitle)

    otherwise
        disp("Option not available");
end










