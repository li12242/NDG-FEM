function ulimit = MomentBDF(mesh, u)
% moment limiter function of Biswas
% 
% Reference: 
% [1]: [Biswas_1994] Parallel, Adaptive Finite Element Methods for Conservation Laws

% moment polynomial coefficient
uh = mesh.Shape.VandMatrix\u;
nOrder = mesh.Shape.nOrder;
eps0=1.0e-8;

ulimit = zeros(size(uh)); 
ulimit(1, :) = uh(1, :);

flag = true(1, mesh.nElement);
for k = (nOrder:-1:1) + 1
    order = k - 1;
    vp = [uh(k-1,2:end), uh(k-1,end)];
    vm = [uh(k-1,1), uh(k-1,1:end-1)];
    norm = LegendreNorm(order);
    norm1 = LegendreNorm(order - 1);
%     norm = 1; norm1 = 1;
    
    ulimit(k, flag) = norm/(2*order - 1)*...
        Utilities.Limiter.Limiter1D.MinmodFun(...
        [ uh(k,flag)*(2*order-1)/norm; ...
        (vp(flag) - uh(k-1,flag) )/norm1; ...
        (uh(k-1,flag) - vm(flag) )/norm1 ]);
    
    ids = (abs(ulimit(k, :) - uh(k, :))<eps0); % the smooth element
    flag(ids) = false;
end% for

ulimit = mesh.Shape.VandMatrix * ulimit;
end% func


function norm = LegendreNorm(k)
norm = sqrt(2/(2*k+1));
end% func