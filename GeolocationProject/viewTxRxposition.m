function viewTxRxposition (bPos, pPos)


num = length(bPos);

if length(pPos )<3
    pPos = [pPos ; 0];   
    bPos = [bPos; zeros(1,size(bPos,2))];
end

figure;
for n = 1: num
    plot(bPos(1,n), bPos(2,n), 'b*','LineWidth',1,'MarkerSize',10); grid on; hold on
end 
plot(pPos(1), pPos(2), 'bo','LineWidth',1,'MarkerSize',10);

figure;
for n = 1: num
    plot3(bPos(1,n), bPos(2,n), bPos(3,n), 'r*','LineWidth',1,'MarkerSize',10); grid on; hold on
end 
plot3(pPos(1), pPos(2), pPos(3), 'bo','LineWidth',1,'MarkerSize',10);