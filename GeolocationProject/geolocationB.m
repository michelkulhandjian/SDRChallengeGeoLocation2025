function [pos] = geolocationB (P, tdoa, M, c, method_flag, p_T_0)

if ~exist('p_T_0','var')
    p_T_0 = 0;
end

if method_flag == 0
    %%% CTLS
    p_1 = P(:,1);
    dummy = P(:,2:M)';
    %d = sqrt( (dummy(:,1)-p_1(1)).^2+(dummy(:,2)-p_1(2)).^2+ (dummy(:,3)-p_1(3)).^2);
    %A = 2*[(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2)), (dummy(:,3)-p_1(3)), c*tdoa];
    if length(p_1) < 3
        diff  = [(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2))];
    else
        diff  = [(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2)), (dummy(:,3)-p_1(3))];
    end
    d = sqrt(sum(diff.^2, 2));
    %d = [c*tdoa
    A = 2*[diff, c*tdoa];
    b = -(c*tdoa).^2 +  sum(diff.^2, 2);
    D2 = diag(d.^(-2));
    ADA = A.'*D2*A;
    cond = (rcond(ADA.'*ADA) <= 1e-12);
    x_lin = pinv(ADA)*A.'*D2*b;
    %x_lin = pinv(A)*b;
    %rmse = norm(p_T-x_lin(1:3))^2;
    if length(p_1) < 3
        pos = x_lin(1:2);
    else
        pos = x_lin(1:3);
    end
elseif method_flag == 1
    %%% OSNM-CTLS Solution
    emin = 0.001; nMax = 1000; ep = 100; n = 1;
    p_1 = P(:,1);
    dummy = P(:,2:M)';
    %d = sqrt( (dummy(:,1)-p_1(1)).^2+(dummy(:,2)-p_1(2)).^2+ (dummy(:,3)-p_1(3)).^2);
    %A = 2*[(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2)), (dummy(:,3)-p_1(3)), c*tdoa];
    if length(p_1) < 3
        diff  = [(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2))];
    else
        diff  = [(dummy(:,1)-p_1(1)), (dummy(:,2)-p_1(2)), (dummy(:,3)-p_1(3))];
    end
    %d = sqrt(sum(diff.^2, 2));
    d = [c*tdoa];
    A = [diff, c*tdoa];
    b = 0.5*(-(c*tdoa).^2 +  sum(diff.^2, 2));
    D2 = diag(d.^(-2));
    ADA = A.'*D2*A;
    cond = (rcond(ADA.'*ADA) <= 1e-12);
    x_lin = pinv(ADA)*A.'*D2*b;
    if length(p_1) < 3
        p = x_lin(1:2);
    else
        p = x_lin(1:3);
    end
    d1 = norm( diff(:,1) );
    while (ep > emin && n <nMax )
        n  = n+1;
        diffpp_1 = p-p_1;
        d = [c*tdoa+ norm(diffpp_1)];
        if length(p_1) < 3
            Jp = [eye(2) (1/d1)*diffpp_1];
        else
            Jp = [eye(3) (1/d1)*diffpp_1];
        end
        JpADA = Jp*ADA*Jp.';
        JpAD  = Jp*A.'*D2*(A*[diffpp_1;norm(diffpp_1)]-b);
        pnew = p- pinv(JpADA)*JpAD;
        epnew = norm(pnew - p)/norm(p);
        if n >2 && epnew > ep && 0
            break;
        else
            ep = epnew;
        end
        p = pnew;
    end
    if length(p_1) < 3
        pos = p(1:2);
    else
        pos = p(1:3);
    end
elseif method_flag == 2 % Taylor Series

    [L, N] = size(P);
    emin = 0.001; nMax = 1000; ep = 100; n = 1;
    if 0
        in_est_error = 20;
        %%% Taylor Series Expansion Solution
        if L < 3
            p_T_0 =  in_est_error*randn(2,1);    %initial estimate with some error (penalty term)
        else
            p_T_0 =  in_est_error*randn(3,1);
        end
    end
    p_T_0 = mean(P,2);
    d = c*tdoa;

    f = zeros(M-1,1);
    del_f = zeros(M-1,length(p_T_0));
    while (ep > emin && n <nMax )
        n  = n+1;

        for ii=2:M
            f(ii-1)=norm(p_T_0-P(:,ii))-norm(p_T_0-P(:,1));
            del_f(ii-1,1) = (p_T_0(1)-P(1,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(1)-P(1,1))*norm(p_T_0-P(:,1))^-1;
            del_f(ii-1,2) = (p_T_0(2)-P(2,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(2)-P(2,1))*norm(p_T_0-P(:,1))^-1;
            if length(p_T_0) >= 3
                del_f(ii-1,3) = (p_T_0(3)-P(3,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(3)-P(3,1))*norm(p_T_0-P(:,1))^-1;
            end
        end
        deltaP = pinv(del_f)*(d-f);
         ep = norm(deltaP);
        p_T_0 = deltaP +p_T_0;

    end
    x_nonlin = deltaP +p_T_0;
    %rmse(k) = norm(p_T-x_nonlin)^2;
    if length(p_T_0) < 3
        pos = x_nonlin(1:2);
    else
        pos = x_nonlin(1:3);
    end
elseif method_flag == 3 % Chan's Ho

    Rnode = P;
    r_i1  = c*tdoa;

    [L, N] = size(Rnode);
    G = []; Y = [];
    for n = 1 : N
        if L <3
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 ;
        else
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 + Rnode(3, n)^2;
        end
    end

    for n = 1 : N
        if n<N
            if L <3
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), r_i1(n)] ];
            else
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), Rnode(3, n+1) - Rnode(3, 1), r_i1(n)] ];
            end
            Y= [Y; k(n+1)-k(1)-r_i1(n)^2];
        end
    end
    G = 2*G;
    %Z_LS=(G'*G)\G'*Y;       % LS Result of Position
    Z_LS=pinv(G'*G)*G'*Y; 
    pos = Z_LS(1:end-1);

elseif method_flag == 4 % Improved Chan's Ho
    X = 0;
    num = 9;  % number of iteration
    emin = 0.001; nMax = 1000; ep = 100; n = 1;

    Rnode = P;
    r_i1  = c*tdoa;

    [L, N] = size(Rnode);
    G = []; Y = [];
    for n = 1 : N
        if L <3
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 ;
        else
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 + Rnode(3, n)^2;
        end
    end

    for n = 1 : N
        if n<N
            if L <3
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), r_i1(n)] ];
            else
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), Rnode(3, n+1) - Rnode(3, 1), r_i1(n)] ];
            end
            Y= [Y; k(n+1)-k(1)-r_i1(n)^2];
        end
    end
    G = 2*G;
    if L <3
        % r1 = Z_LS(3);
        r1 = 1000;
        R11 = [r1];
        XX = G(:,1:2)/2;
        RR = G(:,3)/2;
        RK = -Y;
        % Iteration
        %for nn = 1 : num
        pX = zeros(L,1);
        while (ep > emin && n <nMax )
           n  = n+1;

            [X] = -pinv(XX)*(RR*r1+ (1/2)*RK);
            r1 = norm(X-P(:,1));
            %r11 = sqrt(P(1,1)^2+ P(2,1)^2-2*P(1,1)*X(1)-2*P(2,1)*X(2)+ X(1)^2+X(2)^2);
            %R11 = [R11 r1];
            ep = norm(pX - X);
            pX = X;

        end
        % R11
    end
    pos = X; 
elseif method_flag == 5 % 2-Step Weighting LS solution

    Rnode = P;
    r_i1  = c*tdoa;

    [L, N] = size(Rnode);
    G = []; Y = [];
    for n = 1 : N
        if L <3
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 ;
        else
            k(n) = Rnode(1, n)^2+Rnode(2, n)^2 + Rnode(3, n)^2;
        end
    end

    for n = 1 : N
        if n<N
            if L <3
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), r_i1(n)] ];
            else
                G = [G ; [Rnode(1, n+1) - Rnode(1, 1), Rnode(2, n+1) - Rnode(2, 1), Rnode(3, n+1) - Rnode(3, 1), r_i1(n)] ];
            end
            Y= [Y; k(n+1)-k(1)-r_i1(n)^2];
        end
    end
    G = 2*G;
    SNR = 1;
    % 2 Step ML
    SNR_L=10^(SNR/10);                                              % Linearize SNR

    rr  = 0.5*ones(L)+diag(0.5*ones(1,L));
    cove=4./SNR_L.*rr;
    % Error Covariance Matrix
    %Z=(G'/cove*G)\G'/cove*Y;   % 1 Step Result Z
    invCove = pinv(cove);
    covZ=pinv(G'*invCove*G); 
    Z=covZ*G'*invCove*Y;       % Z Covariance Matirx
    % D Matrix
    if L <3
        D=diag([Z(1)-Rnode(1,1),Z(2)-Rnode(2,1),Z(3)]);
    else
        D=diag([Z(1)-Rnode(1,1),Z(2)-Rnode(2,1),Z(3)-Rnode(3,1), Z(4)]);
    end

    cove1=4*D*covZ*D;                                               % 2 Step Error Covariance Matrix cove1

    if L <3
        G1=[1 0;0 1;1 1];
    else
        G1=[1 0 0;0 1 0; 0 0 1;1 1 1]; % G1 Matrix
    end
    if L <3
        Y1=[(Z(1)-Rnode(1,1))^2;(Z(2)-Rnode(2,1))^2;(Z(3))^2];              % Y1 Matrix
    else
        Y1=[(Z(1)-Rnode(1,1))^2;(Z(2)-Rnode(2,1))^2;(Z(3)-Rnode(3,1))^2; Z(4)^2];              % Y1 Matrix
    end

    % Z2 = (G1'/cove1*G1)\G1'/cove1; 
    invCove1 = pinv(cove1);
    covZ1=pinv(G1'*invCove1*G1);
    Z2 =covZ1*G1'*invCove1; 
    cond = rcond(Z2'*Z2);

    Z_ML=Z2*Y1;                               % 2 Step ML Estimation Result
    Zp=real(sqrt(Z_ML));                      % Extract Position from Result
    if L <3
        Zp=[sign(Z(1)-Rnode(1,1)).*Zp(1);sign(Z(2)-Rnode(2,1)).*Zp(2)]+Rnode(1:2,1);
    else
        Zp=[sign(Z(1)-Rnode(1,1)).*Zp(1);sign(Z(2)-Rnode(2,1)).*Zp(2);sign(Z(3)-Rnode(3,1)).*Zp(3)]+Rnode(:,1);
    end
    pos = Zp;
end
if sum(isnan(pos)) || sum(abs(pos))>1e7
    pos = mean(P,2);
end
end