function phys = Convection1d_Solver(phys)

[rk4a, rk4b, rk4c] = getRK_coefficient;

% set time interval
xm = abs(phys.mesh.x(2) - phys.mesh.x(1)); 
CFL = 0.3;
dt = CFL*min(xm./phys.u(1), xm.^2/sqrt(phys.Dx));

time = 0; outStep = 0;

resT      = zeros(size(phys.var));
FinalTime = phys.ftime;
var       = phys.var;
Np        = phys.mesh.Shape.nNode;
Ne        = phys.mesh.nElement;
while(time < FinalTime)
    
    if time + dt > FinalTime
        dt = FinalTime - time;
    end
    time = time+dt;
    
%     fprintf('Processing: %f ...\n', time./FinalTime)
    
    for INTRK = 1:5
        [rhsT] = Convection1d_RHS(phys, var);
                
        resT = rk4a(INTRK)*resT + dt*rhsT;
        var = var + rk4b(INTRK)*resT;
        
        
%         var = Utilities.Limiter.Limiter1D.BJ1D(phys.mesh,var);
%         [flag, I] = Utilities.Limiter.Limiter1D.DisDetector(mesh, var, u);
%         ind = find(flag);
%         var(:, ind) = temp(:, ind);
%         var = Utilities.Limiter.Limiter1D.MomentBDF(phys.mesh, var);
%         var = Utilities.Limiter.Limiter1D.TVB1D(phys.mesh, var, 1000);
%         var = Utilities.Limiter.Limiter1D.GKTVB1D(phys.mesh, var);
    end% for

%     phys.file.putVarPart('var', [0, 0, outStep],[Np, Ne, 1], var);
%     phys.file.putVarPart('time', outStep, 1, time);
%     outStep = outStep + 1;
%     plot(phys.mesh.x, var); drawnow;
end% while

phys.var = var;
end% func

function [rk4a, rk4b, rk4c] = getRK_coefficient
% rk4c = [0, 1/2, 1/2, 1];
% rk4a = [1/2, 1/2, 1, 1];
% rk4b = [1/6, 1/3, 1/3, 1/6];

%% LSERK
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];
end