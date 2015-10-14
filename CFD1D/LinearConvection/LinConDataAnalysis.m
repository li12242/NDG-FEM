function LinConDataAnalysis
refineOrder
end% func

function refineOrder
% 相同网格，基函数最高基逐渐增大，统计误差
minOrder = 2;
maxOrder = 5;
Order = minOrder:maxOrder; n = numel(Order);
Err = zeros(n, 3); K = 100;
for iorder = 1:n
    [~, Err(iorder,:)] = LinConSetUp(iorder, K);
end
errStr = {'L1', 'L2', 'Linf'};
for ierr =1:3
    subplot(3,1,ierr)
    plot(Order, (Err(:,ierr)), 'ro-')
    xlabel(errStr{ierr})
end
end% func

function refineMesh
% 保持基函数最高阶数不变，逐渐加密网格，统计收敛精度
n = 5; % Number of mesh
K0 = 25; % base element num
nOrder =2; % highest order of base
Err = zeros(n, 3);
K = zeros(n,1);
for imesh =1:n
    K(imesh) = K0*2^(imesh-1);
    [~, Err(imesh,:)] = LinConSetUp(nOrder, K(imesh));
end
figure
for ierr =1:3
    subplot(3,1,ierr)
    plot(log10(K), log10(Err(:,ierr)), 'ro-')
end

% rate: order of accuracy
%     mesh   | L1, L2, L3 |
%    K1-K2   |            |
%    K2-K3   |            |
rate = zeros(n-1, 3);
for ierr = 1:3
    erPlus = Err(2:end, ierr);
    erMius = Err(1:end-1, ierr);
    rate(:,ierr) = log10( erMius./erPlus )./log10(2);
end
% fprintf('\tL1 \t L2 \t Linf \n')
rate
figure
bar(rate, 1.2)
legend('L1', 'L2', 'Linf')
end% func