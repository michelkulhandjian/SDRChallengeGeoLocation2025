function [x1,x2] = quadraRoot(a,b,c)
% Return real root of a*x^2 + b*x + c = 0
x1 = -real((b + (b^2 - 4*a*c)^(1/2))/(2*a));
x2 = -real((b - (b^2 - 4*a*c)^(1/2))/(2*a));

if 0
    r1 = -((b + (b^2 - 4*a*c)^(1/2))/(2*a));
    r2 = -((b - (b^2 - 4*a*c)^(1/2))/(2*a));
    disc = (b^2-4*a*c);
    if r1> 0 && r2>0 % r1<0 && abs(r1) > abs(r2) % disc<= 0
        [disc, r1 , r2]
    end
end
end