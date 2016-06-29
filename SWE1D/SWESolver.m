function [h, q] = SWESolver(phys, ncfile)
% get variables from phys
mesh  = phys.mesh;
bot   = phys.bot;
q     = phys.q; 
h     = phys.h;
ftime = phys.ftime;
dx    = phys.dx;

% Init parameters
time  = 0;
CFL   = 0.2; 
isk   = 0;   % output step
resQ  = zeros(size(q)); 
resH  = zeros(size(h));
% RK coefficient
[rk4a, rk4b, rk4c] = RKcoeff;

% outer time step loop
while(time<ftime)
    % Runge-Kutta residual storage  
    resQ(:) = 0; resH(:) = 0;
    % Estimate time step
    S  = SWESpeed(phys, h, q);
    dt = CFL/S*dx;
    
    if time + dt > ftime
        dt = ftime - time;
    end% if

    fprintf('Processing: %f, dt: %f, wave speed: %f\n',...
        time./ftime, dt, S)

    for INTRK = 1:5
        isWet  = WetDryJudge(phys, mesh, h);
        
        timelocal = time + dt*rk4c(INTRK);       
        [rhsH, rhsQ] = SWERHS(phys, mesh, h, q, bot, isWet);
        
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        resH = rk4a(INTRK)*resH + dt*rhsH;
        
        q = q + rk4b(INTRK)*resQ;
        h = h + rk4b(INTRK)*resH;
        
        [h, q] = PositivePreserving(mesh, h, q);
        
    end
    time = time + dt;
    ncfile.putVarPart('time', isk, 1, time);
    ncfile.putVarPart('h',  [0,0,isk],[mesh.Shape.nNode,mesh.nElement,1], h);
    ncfile.putVarPart('q',  [0,0,isk],[mesh.Shape.nNode,mesh.nElement,1], q);
    
    isk = isk + 1;
end

end% func

function [h, q] = PositivePreserving2(mesh, h, q, bedElva)
h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh, h); 
q = Utilities.Limiter.Limiter1D.MinmodLinear(mesh, q);
q(h<=10^-3) = 0;
h(h<0) = 0;
end% func

function S = SWESpeed(phys, h, q)
% max wave speed
g           = phys.gra;
wetEleFlag  = WetDryJudge(phys, phys.mesh, h);
u           = (q./h) + sqrt(g*h);
S           = max( max( u(:, wetEleFlag)) );
end% func

function [rk4a, rk4b, rk4c] = RKcoeff
% get Runge-Kutta coefficients
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