function [pos] = geolocationA (P, tdoa, M, c, method_flag)

if method_flag == 0
    %%% Linear Solution LS 1
    p_1 = P(:,1);
    dummy = P(:,2:M)';
    if length(p_1) < 3
        A = 2*[(p_1(1)-dummy(:,1)), (p_1(2)-dummy(:,2)), -c*tdoa];
    else
        A = 2*[(p_1(1)-dummy(:,1)), (p_1(2)-dummy(:,2)), (p_1(3)-dummy(:,3)), -c*tdoa];
    end
    b = (c*tdoa).^2 + norm(p_1)^2 - sum((dummy.^2),2); %original one
    %b = [-(dummy(:,1)-p_1(1,1)).^2-(dummy(:,2)-p_1(2,1)).^2-(dummy(:,3)-p_1(3,1)).^2+(c*tdoa).^2];
    x_lin = pinv(A)*b;
    %rmse = norm(p_T-x_lin(1:3))^2;
    if length(p_1) < 3
        pos = x_lin(1:2);
    else
        pos = x_lin(1:3);
    end
elseif method_flag == 1
    %%% Linear Solution LS 2
    p_1 = P(:,1);
    dummy = P(:,2:M)';
    if length(p_1) < 3
        A = -[(p_1(1)-dummy(:,1)), (p_1(2)-dummy(:,2)), -c*tdoa];
    else
        A = -[(p_1(1)-dummy(:,1)), (p_1(2)-dummy(:,2)), (p_1(3)-dummy(:,3)), -c*tdoa];
    end
    %b = (c*tdoa).^2 + norm(p_1)^2 - sum((dummy.^2),2); %original one
    if length(p_1) < 3
        b = .5*[(dummy(:,1)-p_1(1,1)).^2+(dummy(:,2)-p_1(2,1)).^2-(c*tdoa).^2];
    else
        b = .5*[(dummy(:,1)-p_1(1,1)).^2+(dummy(:,2)-p_1(2,1)).^2+(dummy(:,3)-p_1(3,1)).^2-(c*tdoa).^2];
    end
    x_lin = pinv(A)*b;
    if length(p_1) < 3
        pos = x_lin(1:2);
    else
        pos = x_lin(1:3);
    end

elseif method_flag == 2  % good case
    root_sel = 1;
    [L, N] = size(P);
    p_1 = P(:,1);
    W = diag(ones(1,N-1));
    PP = P(:,2:M)-p_1;
    S = PP';
    R = (sum(PP.^2))';
    delta = (c*tdoa);  % d vector
    RR = R-delta.^2;  % deltaV
    Sws = pinv(S'*W*S)*S'*W;
    Ps = S*Sws;
    PsOrth = eye(N-1)-Ps;
    if 0
        % Option 1
        Rs = abs(delta'*PsOrth*PsOrth*RR/(2*delta'*PsOrth*PsOrth*delta));
    else
        % Option 2
        % a = 4 - 4*(delta'*Sws'*Sws*delta); b = 4*(delta'*Sws'*Sws*RR);
        % c1 = -RR'*Sws'*Sws*RR;
        a1 = delta'*Sws'*Sws*delta;
        a = 4 - 4*(a1); b = 4*(delta'*Sws'*Sws*RR);
        c1 = -RR'*Sws'*Sws*RR; d1 = b^2-4*a*c1;
        [r_1,r_2] = quadraRoot(a,b,c1);

        if a1 < 1  % there exist unique real positive root - r_2
            Rs = r_2;
        elseif a1 > 1
            if d1 >= 0
                if b < 0
                    Rs = r_1;  % there exist unique real positive root - r_1
                else % b >= 0
                    % it is bad geometry and we should filter out.
                    % The resulting r_1 and r_2 either both belong to true
                    % hyperbola intersection or dont belong to solution.
                    % Also, r_1 and r_2 to the same solution. 
                    % There is an ambiguity In this case, at least 4 sensors
                    % in 2D and 5 sensors in 3D are needed to resolve the
                    % ambiguity.
                    % The are some possibilites depending on the geometry,
                    % assumption of Tx location and channel conditions:
                    % a. No estimate location by assigning nan vector.
                    % b. Select the shortest positive root among r_1 and r_2.
                    % c. We can check r_1 if positive the set Rs = r_1 if not Rs = r_2;

                    %% a.
                    if root_sel == 0
                        Rs = nan;
                    end
                    %% b.
                    if root_sel ==1
                        if (r_1 < r_2)
                            if (r_1 > 0)
                                Rs = r_1;
                            else
                                Rs = r_2;
                            end
                        else
                            if (r_2 > 0)
                                Rs = r_2;
                            else
                                Rs = r_1;
                            end
                        end
                    end
                    %% c.
                    if 0
                        if r_1 > 0
                            Rs = r_1;
                        else
                            Rs = r_2;
                        end
                    end
                end
            else % d1 < 0
                % Results in imaginary roots due to bad geometry we should filter out 
                % However, the real parts of both r_1 and r_2 are the same
                % and belong to true solution.

                % The are some possibilites depending on the geometry,
                % assumption of Tx location and channel conditions:
                % a. No estimate location by assigning nan vector.
                % b. Select the r_1 since both roots are equal.
                %% a.
                if root_sel == 0
                    Rs = nan;
                end
                %% b.
                if root_sel == 1
                    Rs = r_1;
                end
            end
        end % if a1 < 1

    end % if 0 or the other

    % pos1 = 0.5*Sws*(RR-2*r_1*delta);
    % pos2 = 0.5*Sws*(RR-2*r_2*delta);
    % erV1 = RR-2*r_1*delta -2*S*pos1;
    % erV2 = RR-2*r_2*delta -2*S*pos2;
    % if norm(erV1)< norm(erV2)
    %     Rs = r_1;
    % else
    %     Rs = r_2;
    % end


    % if r_1 > 0 && r_2 > 0
    %     if r_1 < r_2
    %         Rs = r_1;
    %     else
    %         Rs = r_2;
    %     end
    % elseif  r_1 > 0
    %     Rs = r_1;
    % else
    %     Rs = r_2;
    % end
    pos = p_1 + 0.5*Sws*(RR-2*Rs*delta);

elseif method_flag == 3
    [L, N] = size(P);
    p_1 = P(:,1);
    W = diag(ones(1,L));
    PP = P(:,2:M)-p_1;
    S = PP';
    R = (sum(PP.^2))';
    delta = (c*tdoa);       % d
    RR = 0.5*(R-delta.^2);  % z
    Sws = pinv(S'*W*S)*S'*W;
    aa = Sws*RR;
    bb = Sws*delta;
    a1 = aa(1); a2 = aa(2);
    b1 = bb(1); b2 = bb(2);
    Rs = sqrt((a1*b1+a2*b2)^2 - (b1^2+b2^2-1)*(a1^2+a2^2))/(b1^2+b2^2-1);
    pos = aa-Rs*bb;

elseif method_flag == 4
    [L, N] = size(P);
    p_1 = P(:,1); x1 = P(1,1); x2 = P(1,2); x3 = P(1,3);
    y1 = P(2,1); y2 = P(2,2); y3 = P(2,3);
    delta = (c*tdoa);
    m1 = 0;
    m2 = 0.5*(norm(P(:,1))^2 - norm(p_1)^2 - delta(1)^2);
    alpha = (delta(1)*(x3-x1) - delta(2)*(x2-x1))/(-delta(1)*(y3-y1)+delta(2)*(y2-y1));
    beta = (delta(2)*m1-delta(1)*m2)/(-delta(1)*(y3-y1)-delta(2)*(y2-y1));

    a = (1+alpha^2) - ((x2-x1)+(y2-y1)*alpha)^2;
    b = -2*(x1+(y1-beta)*alpha-((x2-x1)+(y2-y1)*alpha)*((y2-y1)*beta-m1));
    c = (y1-beta)^2 + x1^2 - (((y2-y1)*beta-m1)/delta(1)^2)^2;

    r1 = (-b + sqrt(b^2-4*a*c))/2/a;
    r2 = (-b - sqrt(b^2-4*a*c))/2/a;
    if r1> 0
        xs = r1;
    else
        xs = r2;
    end
    ys = alpha*xs+ beta;

    pos = [xs ys]';


elseif method_flag == 5
    [L, N] = size(P);
    d1 = 1000;
    p_1 = P(:,1);
    W = diag(ones(1,L));
    PP = P(:,2:M)-p_1;
    S = [zeros(1,L);-PP'];
    %R = (sum(PP.^2))';
    delta = (c*tdoa);
    G = [S, [1;delta]];
    h = [d1;0.5*(norm(p_1)^2-(sum(P(:,2:M).^2))'+delta.^2)];
    W = diag(ones(1,L+1));
    wp = inv(G'*W*G)*G'*W*h;
    W2 = diag(1./([[0, delta']+wp(end)].^2) );
    wp2 = inv(G'*W2*G)*G'*W2*h;
    pos = wp2(1:end-1);

elseif method_flag == 6
    bPos = P;
    xp = mean(bPos(1,:));
    yp = mean(bPos(2,:));

    x0 = bPos(1, 1); y0 = bPos(2, 1);
    x1 = bPos(1, 2); y1 = bPos(2, 2);
    x2 = bPos(1, 3); y2 = bPos(2, 3);
    rangeDif =c*tdoa;

    F = functions(x0, y0, x1, y1, x2, y2, rangeDif(1), rangeDif(2), rangeDif(3));
    options = optimset('Display', 'off');
    [vOpt, ~] = lsqnonlin(F, [xp, yp], [], [], options);
    x = vOpt(1);
    y = vOpt(2);

    pos = [x y]';

end
if 0 
if sum(isnan(pos)) || sum(abs(pos))>1e7
    pos = mean(P,2);
end
end

end
% fails = sum(rmse > fail_thr^2);
% RMSE(m) = sqrt(mean(rmse(rmse < fail_thr^2)));

function fun = functions(x0, y0, x1, y1, x2, y2, d01, d02, d12)

%fun = @(x, y) fn(x, y)
%function A = fn(x, y)
fun = @(xy) fn(xy)
    function A = fn(xy)
        %function [a1,b1,c1] = fn(xy)
        x = xy(1);
        y = xy(2);
        a = sqrt((x - x1)^2 + (y - y1)^2) - sqrt((x - x0)^2 + (y - y0)^2) - d01;
        b = sqrt((x - x2)^2 + (y - y2)^2) - sqrt((x - x0)^2 + (y - y0)^2) - d02;
        c = sqrt((x - x2)^2 + (y - y2)^2) - sqrt((x - x1)^2 + (y - y1)^2) - d12;
        A = [a, b, c];
        [a1,b1,c1] = deal(a, b, c);

    end

end