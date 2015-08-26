function [h, q] = SWESolver(mesh, h, q, bedElva, FinalTime)
% 1D shallow water equation 
% Purpose: time setpping
% REFERENCE:
% [1]: [Xing_2010] Positivity-preserving high order well-balanced
%      discontinuous Galerkin methods for the shallow water equations.
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
     
time = 0;
% Runge-Kutta residual storage  
resQ = zeros(size(q)); resH = zeros(size(h));

% compute time step size
xmin = min(abs(mesh.x(1,:)-mesh.x(2,:)));

CFL=0.20; 
% is_Camera_on = 1;	% 否进行摄像
% if is_Camera_on
%     writerObj = VideoWriter('SWE1D.avi');
%     writerObj.FrameRate=50;
%     open(writerObj);	% 设定好储存的帧速率以后才能打开
% end
% outer time step loop 
while(time<FinalTime)
    lamda = SWESpeed(h, q);
    
    dt = CFL/max(lamda(:))*xmin;
    for INTRK = 1:5
        
        [rhsH, rhsQ] = SWERHS(mesh, h, q, bedElva);
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        resH = rk4a(INTRK)*resH + dt*rhsH;
        q = q + rk4b(INTRK)*resQ;
        h = h + rk4b(INTRK)*resH;
        [h, q] = PositivePreserving(mesh, h, q, bedElva);

        % plot
%         subplot(2,1,1)
%     if time > 20
%         plot(mesh.x, h+bedElva, 'o-'); drawnow
%     end
%         subplot(2,1,2)
%         plot(mesh.x, q, 'o-'); drawnow
        
%         if is_Camera_on
%             frame = getframe(gcf);
%             writeVideo(writerObj,frame);     % writeVideo 语句并不返回任何参数
%         end
    end
    % Increment time
    time = time+dt;
    
end
% if is_Camera_on
%     close(writerObj);
% end
end% func

function [h, q] = PositivePreserving(mesh, h, q, bedElva)
% Positivity-preserving methods
% REFERENCE:
% [1]: Xing(2010)
eta = h + bedElva;
hlim = Utilities.SlopeLimitN(mesh, h); q = Utilities.SlopeLimitN(mesh, q);
etalim = Utilities.SlopeLimitN(mesh, eta);

hmin = min(h); etaLimit = find(hmin>0); hLimit = find(hmin<0);
% for hmin>0, TVB limiter is based on (h+b, q)
h(:, etaLimit) = etalim(:, etaLimit) - bedElva(:, etaLimit);
% for hmin<0, TVB limiter is based on (h, q)
h(:, hLimit) = hlim(:, hLimit);

if(~isempty(hLimit))
    hdelta = 10^-12; % scheme min depth
    % Compute cell averages
    hmean = CellMean(mesh, h(:, hLimit)); hmean(hmean< 10^-12) = 10^-12;
    theta = min(1, hmean./(hmean - hmin(hLimit) + hdelta));
    hmean = repmat(hmean, mesh.Shape.nNode, 1);
    theta = repmat(theta, mesh.Shape.nNode, 1);
    h(:, hLimit) = theta.*(h(:, hLimit) - hmean) + hmean;
end% if

% eliminate the flux in near dry regions
q(h<10e-6) = 0;
end% func

function hmean = CellMean(mesh, h)
Np = mesh.Shape.nNode;
uh = mesh.Shape.VandMatrix\h; uh(2:Np,:)=0;
uavg = mesh.Shape.VandMatrix*uh; hmean = uavg(1,:);
end% func

function lambda = SWESpeed(h, q)
% max wave speed
TOL = 0; g = 9.8;
flag = (h>TOL);
u = (q./h) + sqrt(g*h);
lambda = max(u(flag));
end% func