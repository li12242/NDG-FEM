function var = Convection2DSolver(mesh, var, FinalTime, u, v, outfile)
% 2D convection problem
% RK time stepping

[rk4a, rk4b, rk4c] = RK4Coeff;
time = 0;
% Filt = Fliter(mesh.Shape, mesh.Shape.nOrder, 0.9);

resVar = zeros(size(var));
% compute time step size
xmin = min(sqrt((mesh.x(1,:)-mesh.x(2,:)).^2 + (mesh.y(1,:) - mesh.y(2,:)).^2 ));
CFL=0.30;  
un = max(max( sqrt(u.^2 + v.^2) ));
dt = CFL/un*xmin; outStep = 0;

while(time < FinalTime)
    
    if time + dt > FinalTime
        dt = FinalTime - time;
    end
    time = time+dt;
    
    fprintf('Processing: %f ...\n', time./FinalTime)
    
    for INTRK = 1:5
        rhsVar = Convection2DRHS(mesh, var, time, u, v);
        
%         % filter residual
%         rhsVar = Filt*rhsVar;
        
        resVar = rk4a(INTRK)*resVar + dt*rhsVar;
        var = var + rk4b(INTRK)*resVar;
        
        temp = Utilities.Limiter.Limiter2D.BJ2D(mesh, var);
%         temp = Utilities.Limiter.Limiter2D.JKTA_tri(mesh, var);
        [flag, I] = Utilities.Limiter.Limiter2D.DisDetector(mesh, var, u, v);
        ind = find(flag);
        var(:, ind) = temp(:, ind);
        
    end% for
    outfile.putVarPart('var', [0, outStep], [mesh.nNode, 1], var);
    outfile.putVarPart('time', outStep, 1, time);
    
    outStep = outStep + 1;
%     plot3(mesh.x(mesh.vmapP),mesh.y(mesh.vmapP), var(mesh.vmapM)); drawnow
end% while
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