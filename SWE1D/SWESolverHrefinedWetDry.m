function [h, q] = SWESolverHrefinedWetDry(physics, ncfile)
% time setpping of 1D shallow water equation 
% 

mesh = physics.getVal('mesh');
bedElva = physics.getVal('bedElva');

time = 0;
q = physics.getVal('flux'); 
h = physics.getVal('height');
[rk4a, rk4b, rk4c] = RKcoeff;

% compute time step size
xmin = min(abs(mesh.x(1,:)-mesh.x(2,:)));
CFL=0.3; outstep = 0;
FinalTime = physics.getVal('FinalTime');


% eliminate zero depth in wet cell
isWet = WetDryJudge(mesh, h, physics);
[h, q] = PositivePreserving(mesh, h, q, bedElva, isWet);

% outer time step loop
while(time<FinalTime)
    lamda = SWESpeed(h, q, physics, isWet);
    dt = CFL/lamda*xmin;
    
    if dt < 0.005 
        dt = 0.005;
    end
    
    % Increment time
    if time + dt > FinalTime
        time = FinalTime;
        dt = FinalTime - time;
    else
        time = time + dt;
    end% if
    

    fprintf('Processing: %f, dt: %f, wave speed: %f\n',...
        time./FinalTime, dt, lamda)
    
%     if time/FinalTime > .126160
%         keyboard
%     end% if
%     if outstep > 577
%         keyboard
%     end

    refineflag = RefinedCellIdentify(mesh, h, physics, isWet);
    [new_mesh, h1, q1, new_bedElva, localEleIndex] =...
        Hrefine1D(mesh, refineflag, h, q, bedElva, physics);
    isWet = WetDryJudge(new_mesh, h1, physics);
    % Runge-Kutta residual storage  
    resQ = zeros(size(q1)); resH = zeros(size(h1));

    for INTRK = 1:5
        
%         subplot(3,1,1); 
%         plot(new_mesh.x, h1+new_bedElva, '-b.', new_mesh.x, new_bedElva, 'k'); 
%         subplot(3,1,2); plot(new_mesh.x, q1, '-r'); 
%         u = q1./h1; u(h1<1e-3) = 0;
%         subplot(3,1,3); plot(new_mesh.x, u, '-b.');
%         drawnow;
        
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
    StoreVar(ncfile, h, q, time, lamda, outstep)
    outstep = outstep + 1;
end

end% func

function [h, q] = PositivePreserving2(mesh, h, q, bedElva)
h = Utilities.Limiter.SlopeLimit1(mesh, h); 
q = Utilities.Limiter.SlopeLimit1(mesh, q);
% q(h<=10^-3) = 0;
h(h<0) = 0;
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