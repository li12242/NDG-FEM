function ulimit = SlopeLimitBiswas(mesh, u)
% Reference: 
% [1]: [Biswas_1994] Parallel, Adaptive Finite Element Methods for 
%       Conservation Laws.pdf

uh = mesh.Shape.VandMatrix\u;
nMod = size(uh, 1);

ulimit = zeros(size(uh)); 
ulimit(1, :) = uh(1, :);

for k = nMod:-1:2
    vp = [uh(k-1,2:end), uh(k-1,end)];
    vm = [uh(k-1,1), uh(k-1,1:end-1)];
    ulimit(k, :) = LegendreNorm(k+1)/(2*k+1)*...
        Utilities.Limiter.minmod([ uh(k,:)*(2*k+1)/LegendreNorm(k+1); ...
        (vp - uh(k-1,:) )/LegendreNorm(k); ...
        (uh(k-1,:) - vm )/LegendreNorm(k) ]);
end% for

ulimit = mesh.Shape.VandMatrix * ulimit;
end% func


function norm = LegendreNorm(k)
norm = sqrt(2/(2*k+1));
end% func