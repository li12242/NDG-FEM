function var = Convection1DSolver(mesh, var, FinalTime, a)


[rk4a, rk4b, rk4c] = getRK_coefficient;

% set time interval
xm = abs(mesh.x(2) - mesh.x(1)); CFL = 0.5;
dt = CFL*(xm./a);

time = 0; outStep = 0;

resVar = zeros(size(var));

while(time < FinalTime)
    
    if time + dt > FinalTime
        dt = FinalTime - time;
    end
    time = time+dt;
    
    fprintf('Processing: %f ...\n', time./FinalTime)
    
    for INTRK = 1:5
        rhsVar = Convection1DRHS(mesh, var, a);
                
        resVar = rk4a(INTRK)*resVar + dt*rhsVar;
        var = var + rk4b(INTRK)*resVar;
        
    end% for
%     StoreVar(outfile, mesh, var, I, time, outStep);
    outStep = outStep + 1;
    
end% while

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