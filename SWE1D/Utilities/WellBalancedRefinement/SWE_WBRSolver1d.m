%% SWE_WBRSolver1d
% Well-balanced refinement solver for one dimensional SWE.
function [h, q] = SWE_WBRSolver1d(physics, ncfile)
% time setpping of 1D shallow water equation 
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

% eliminate zero depth in wet cell
isWet = WetDryJudge(mesh, h, physics);
% [h, q] = PositivePreserving(mesh, h, q, bedElva, isWet);
lamda = SWESpeed(h, q, physics, isWet);
StoreVar(ncfile, h, q, time, lamda, isk)
% outer time step loop
while(time<FinalTime)
    lamda = SWESpeed(h, q, physics, isWet);
    dt = CFL/lamda*xmin;
    
    % Increment time
    if time + dt > FinalTime
        time = FinalTime;
        dt = FinalTime - time;
    else
        time = time + dt;
    end% if
    

    fprintf('Processing: %f, dt: %f\n', time./FinalTime, dt)
    
    refineflag = RefinedCellIdentify(mesh, h, physics, isWet);
    [new_mesh, h1, q1, new_bedElva, localEleIndex] =...
        Hrefine1D(mesh, refineflag, h, q, bedElva, physics);
    isWet = WetDryJudge(new_mesh, h1, physics);
    % Runge-Kutta residual storage  
    resQ = zeros(size(q1)); resH = zeros(size(h1));

    for INTRK = 1:5
%         timelocal = time + dt*rk4c(INTRK);

        [rhsH, rhsQ] = SWERHS(new_mesh, h1, q1, new_bedElva, isWet, physics);
        
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        resH = rk4a(INTRK)*resH + dt*rhsH;
        
        q1 = q1 + rk4b(INTRK)*resQ;
        h1 = h1 + rk4b(INTRK)*resH;
   
        [h1, q1] = PositivePreserving(new_mesh, h1, q1, new_bedElva);
    end
    [h, q] = Hcombine1D(localEleIndex, h1, q1, new_mesh, refineflag);
    isWet = WetDryJudge(mesh, h, physics);
    StoreVar(ncfile, h, q, time, lamda, isk)
    isk = isk + 1;
end

end% func

function lambda = SWESpeed(h, q, physics, isWet)
% max wave speed
g = physics.getVal('gravity');
u = (q./h) + sqrt(g*h);
lambda = max( max(u(:, isWet)) );
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