function var = Convection2d_Solver(phys)
% 2D convection problem

% parameters
mesh    = phys.mesh;
shape   = mesh.Shape;
Ne      = mesh.nElement;
Np      = shape.nNode;
var     = phys.var;
ftime   = phys.ftime;
u       = phys.u;
v       = phys.v;
ncfile  = phys.file;
% time step
CFL     = 0.30;
un      = max(max( sqrt(u.^2 + v.^2) ));
dt      = CFL*phys.dx/un;

% RK time stepping parameters
[rk4a, rk4b, rk4c] = RK4Coeff;

% variables
time    = 0;
resVar  = zeros(size(var));
contour = 0;

% store initial result
ncfile.putVarPart('var', [0, 0, contour],[Np, Ne, 1], var);
ncfile.putVarPart('time', contour, 1, time);
contour = contour + 1;

while(time < ftime)
    
    if time + dt > ftime
        dt = ftime - time;
    end
    time = time+dt;
    fprintf('Processing: %f ...\n', time./ftime)
    
    for INTRK = 1:5
        rhsVar = Convection2d_RHS(mesh, var, u, v);
                
        resVar = rk4a(INTRK)*resVar + dt*rhsVar;
        var    = var + rk4b(INTRK)*resVar;
        
%         ind    = Utilities.Limiter.Limiter2D.KXRCF2d(mesh, var, u, v, 0.05);
        ind    = Utilities.Limiter.Limiter2D.TVB_detector2d(mesh, var, 0.1);
        flag   = ind > 0;
%         figure('Position', [683     1   557   984]); 
%         subplot(2,1,1); hold on;
%         plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), var(mesh.vmapM), 'k-');
%         plot3(mesh.x(:, flag), mesh.y(:, flag), var(:, flag), 'ro');
        
%         varlim = Utilities.Limiter.Limiter2D.VB2d(mesh, var);
%         varlim = Utilities.Limiter.Limiter2D.VB2d(mesh, var);
        varlim = Utilities.Limiter.Limiter2D.TVB2d(mesh, var);
%         varlim = Utilities.Limiter.Limiter2D.JKTA_tri(mesh, var);
%         var = Utilities.Limiter.Limiter2D.JKTA_quad(mesh, var);
%         var = Utilities.Limiter.Limiter2D.BJ2(mesh, var, 1);
%         var = Utilities.Limiter.Limiter2D.BJLoc2(mesh, var, 1);
%         var = Utilities.Limiter.Limiter2D.HWENO2d(mesh, var);
        
        var(:, flag) = varlim(:, flag);
        
%         subplot(2,1,2)
%         plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), var(mesh.vmapM), 'k-');
    end% for
    
    ncfile.putVarPart('var', [0, 0, contour],[Np, Ne, 1], var);
    ncfile.putVarPart('time',contour,1,time);
    
    contour = contour + 1;
end% while
    
ncfile.CloseFile();
end% func

function [rk4a, rk4b, rk4c] = RK4Coeff
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
end% func

function F = Fliter(Shape, Nc, frac)

filterdiag = ones(Shape.nNode, 1);

% build exponential filter
sk = 1; N = Shape.nOrder;
for i=0:N
  for j=0:N-i
    if (i+j>=Nc)
      filterdiag(sk) = frac;
    end
    sk = sk+1;
  end
end

F = ( Shape.VandMatrix*diag(filterdiag) )/Shape.VandMatrix;
end% func