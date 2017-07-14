function ulimit = MomentBDF(mesh, u)
% moment limiter function of Biswas
% 
% Reference: 
% 1. Biswas(1994)

% moment polynomial coefficient
uh = mesh.Shape.VandMatrix\u;
nOrder = mesh.Shape.nOrder;
eps0 = 1.0e-8;

ulimit = uh; % limit value initialization

flag = true(1, mesh.nElement);
for k = (nOrder:-1:1) + 1
    order = k - 1; v = uh(k-1,:);
    
    vm = v(mesh.EToE(:,1)');
    vp = v(mesh.EToE(:,2)');
   
    norm = LegendreNorm(order)*LegendreNorm(order-1);
    
    ulimit(k, flag) = 1/norm*...
        Utilities.Limiter.Limiter1D.MinmodFun(...
        [ uh(k, flag)*norm; ...
        (vp(flag) - v(flag) ); ...
        (v(flag) - vm(flag) ) ]);
        
    ids = (abs(ulimit(k, :) - uh(k, :))<eps0); % the smooth element
    flag(ids) = false;
end% for

ulimit = mesh.Shape.VandMatrix * ulimit;

end% func


function norm = LegendreNorm(k)
norm = sqrt(2/(2*k+1));
end% func