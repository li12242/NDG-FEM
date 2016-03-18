function [h, q] = SWESolver(physics, ncfile)
% time setpping of 1D shallow water equation 
% 
% REFERENCE:
% [1]: [Xing_2010] Positivity-preserving high order well-balanced
%      discontinuous Galerkin methods for the shallow water equations.

mesh = physics.getVal('mesh');
bedElva = physics.getVal('bedElva');

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
     
time = 0;
q = physics.getVal('flux'); h = physics.getVal('height');
% Runge-Kutta residual storage  
resQ = zeros(size(q)); resH = zeros(size(h));

% compute time step size
xmin = min(abs(mesh.x(1,:)-mesh.x(2,:)));
CFL=0.3;
FinalTime = physics.getVal('FinalTime');
outstep = 0;
lamda = SWESpeed(h, q);
dt = CFL/lamda*xmin;

% outer time step loop
while(time<FinalTime)
    lamda = SWESpeed(h, q);
%     dt = CFL/max(lamda(:))*xmin;
    % Increment time
    if time + dt > FinalTime
        time = FinalTime;
        dt = FinalTime - time;
    else
        time = time + dt;
    end% if
    
    if lamda*dt > CFL*xmin
        dt = dt/2;
    elseif lamda*dt < CFL*xmin/4
        dt = dt*2;
    end%if
%     if time >= 26.6555 - 3*dt
%         keyboard
%     end
    
    fprintf('Processing: %f, dt: %f, wave speed: %f\n',...
        time./FinalTime, dt, lamda)
    for INTRK = 1:5
%         try
        timelocal = time + dt*rk4c(INTRK);
        [rhsH, rhsQ] = SWERHS(mesh, h, q, bedElva);
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        resH = rk4a(INTRK)*resH + dt*rhsH;
        
        q = q + rk4b(INTRK)*resQ;
        h = h + rk4b(INTRK)*resH;

%         plot(mesh.x, h+bedElva, '-b.', mesh.x, bedElva, 'k');
%         drawnow;
    end
    StoreVar(ncfile, h, q, time, lamda, outstep)
    outstep = outstep + 1;
end
% if is_Camera_on
%     close(writerObj);
% end
end% func

function [h, q] = PositivePreserving2(mesh, h, q, bedElva)
h = Utilities.Limiter.SlopeLimit1(mesh, h); 
q = Utilities.Limiter.SlopeLimit1(mesh, q);
q(h<=10^-3) = 0;
h(h<0) = 0;
end% func

function [h, q] = PositivePreserving(mesh, h, q, bedElva)
% Positivity-preserving methods
% REFERENCE:
% [1]: Xing(2010)
hDelta = 0.000; % scheme min depth
limitH = 0.001;
eta = h + bedElva;
q = Utilities.Limiter.SlopeLimitN(mesh, q);

% 1. limiter on eta
etalim = Utilities.Limiter.SlopeLimitN(mesh, eta);
% 2. limiter on h
hlim = Utilities.Limiter.SlopeLimitN(mesh, h); 
% for hmin<0, TVB limiter is based on (h, q)
h = hlim;
% for hmin>0, TVB limiter is based on (h+b, q)
hlim = etalim - bedElva;
hmin = min(hlim); etaLimit = find(hmin > limitH);
h(:, etaLimit) = etalim(:, etaLimit) - bedElva(:, etaLimit);

hmin = min(h); hLimit = find(hmin<hDelta);
if(~isempty(hLimit))
%     hdelta = 10^-12; % scheme min depth
    % Compute cell averages
    hmean = CellMean(mesh, h(:, hLimit)); hmean(hmean<= hDelta) = hDelta;
    theta = min(1, hmean./(hmean - hmin(hLimit) + hDelta));
    hmean = repmat(hmean, mesh.Shape.nNode, 1);
    theta = repmat(theta, mesh.Shape.nNode, 1);
    h(:, hLimit) = theta.*(h(:, hLimit) - hmean) + hmean;
end% if
% temp
% h(h<hDelta) = hDelta;
% eliminate the flux in near dry regions
q(h<=0) = 0;
% q(h<2*10^-2) = 0;
end% func

function hmean = CellMean(mesh, h)
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func

function lambda = SWESpeed(h, q)
% max wave speed
TOL = 0.03; g = 9.8;
flag = (h>TOL);
u = (q./h) + sqrt(g*h);
lambda = max(u(flag));
end% func