function [t]=calc_Delay(Rnode, target)
% Get the speed of light in meters per second
c = physconst('LightSpeed');

N = length(Rnode);

% Calculate Actual Delays To Each Anchor Node
for n = 1 : N
     t(n)= norm(Rnode(:,n) -target)/c;
end
