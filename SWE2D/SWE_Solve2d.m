function phys = SWE_Solve2d(phys, outfile)
%% Initializion of variables
h         = phys.h;    % water depth
qx        = phys.qx;   % water flux
qy        = phys.qy;   % water flux
bot       = phys.bot;  % topography
mesh      = phys.mesh;
FinalTime = phys.ftime;
contour   = 0;
outStep   = 1;
stride    = 1;  % output stride
dx        = phys.dx;   % mesh length
dtm       = phys.dt;   % ideal time step
% Parameters
Np        = mesh.Shape.nNode; % No. of nodes in each element
Ne        = mesh.nElement; % No. of element
% Filt      = Fliter(mesh.Shape, mesh.Shape.nOrder, 0.9);
CFL       = 0.3;

% 5-stage Runge-Kutta coefficients
[rk4a, rk4b, rk4c] = RK45coef;
time = 0;
% Runge-Kutta residual storage  
resH  = zeros(size(h));
resQx = zeros(size(qx));
resQy = zeros(size(qy));

%% RK time stepping
tic
while(time<FinalTime)
    s  = SWE_PredictWaveSpeed2d(phys, h, qx, qy);
    dt = CFL*dx/s;
    
    if(dt<dtm)
        dt = dtm;
    end% if
    if(time+dt>FinalTime)
        dt = FinalTime - time;
    end% if
    
    % zeros Runge-Kutta residual storage
    resQx(:) = 0; resQy(:) = 0; resH(:) = 0;
    for INTRK = 1:5
        timeloc = time + rk4c(INTRK)*dt;
        [rhsH, rhsQx, rhsQy] = SWE_RHS2d(phys, mesh, h, qx, qy, timeloc);
        
        % filter
%         rhsH  = Filt*rhsH;
%         rhsQx = Filt*rhsQx;
%         rhsQy = Filt*rhsQy;
        
        resH  = rk4a(INTRK)*resH  + dt*rhsH;
        resQx = rk4a(INTRK)*resQx + dt*rhsQx;
        resQy = rk4a(INTRK)*resQy + dt*rhsQy;
        h     = h  + rk4b(INTRK)*resH;
        qx    = qx + rk4b(INTRK)*resQx;
        qy    = qy + rk4b(INTRK)*resQy;
        
        %% Slope limiter
        % In ParabolicBowl test case, SL2 works well for quadrilateral 
        % mesh, while SLLoc2 works better for triangle mesh.
        
%         zeta = h+bot;
%         zeta = Utilities.Limiter.Limiter2D.BJ2(mesh, zeta, 1);
%         h    = zeta-bot;
%         h  = Utilities.Limiter.Limiter2D.VB2d_VA(mesh, h);
%         qx = Utilities.Limiter.Limiter2D.VB2d_VA(mesh, qx);
%         qy = Utilities.Limiter.Limiter2D.VB2d_VA(mesh, qy);

%         h  = Utilities.Limiter.Limiter2D.VB2d_JK(mesh, h);
%         qx = Utilities.Limiter.Limiter2D.VB2d_JK(mesh, qx);
%         qy = Utilities.Limiter.Limiter2D.VB2d_JK(mesh, qy);

%         h  = Utilities.Limiter.Limiter2D.VB2d_HWENO(mesh, h);
%         qx = Utilities.Limiter.Limiter2D.VB2d_HWENO(mesh, qx);
%         qy = Utilities.Limiter.Limiter2D.VB2d_HWENO(mesh, qy);

%         h  = Utilities.Limiter.Limiter2D.TVB_tri2d(mesh, h,  0.);
%         qx = Utilities.Limiter.Limiter2D.TVB_tri2d(mesh, qx, 0.);
%         qy = Utilities.Limiter.Limiter2D.TVB_tri2d(mesh, qy, 0.);
        
%         h  = Utilities.Limiter.Limiter2D.TVB_quad2d(mesh, h,  2);
%         qx = Utilities.Limiter.Limiter2D.TVB_quad2d(mesh, qx, 2);
%         qy = Utilities.Limiter.Limiter2D.TVB_quad2d(mesh, qy, 2);

        h  = Utilities.Limiter.Limiter2D.BJ2(mesh, h,  1);
        qx = Utilities.Limiter.Limiter2D.BJ2(mesh, qx, 1);
        qy = Utilities.Limiter.Limiter2D.BJ2(mesh, qy, 1);

%         h  = Utilities.Limiter.Limiter2D.JKTA_tri(mesh, h);
%         qx = Utilities.Limiter.Limiter2D.JKTA_tri(mesh, qx);
%         qy = Utilities.Limiter.Limiter2D.JKTA_tri(mesh, qy);
        
%         h  = Utilities.Limiter.Limiter2D.JKTA_quad(mesh, h);
%         qx = Utilities.Limiter.Limiter2D.JKTA_quad(mesh, qx);
%         qy = Utilities.Limiter.Limiter2D.JKTA_quad(mesh, qy);
        %% Positive-preserving limiter
        [h, qx, qy] = SWE_PositivePreserving2d(phys, mesh, h, qx, qy);
        
%         figure('Position', [627     1   561   984]);
%         subplot(3,1,1);
%         plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), h(mesh.vmapM), 'k.-');
%         view([60, 33])
%         subplot(3,1,2);
%         plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), qx(mesh.vmapM), 'k.-');
%         subplot(3,1,3);
%         plot3(mesh.x(mesh.vmapM), mesh.y(mesh.vmapM), qy(mesh.vmapM), 'k.-');
        
    end
    
    % Increment time
    time = time+dt;
    
    outfile.putVarPart('time', outStep, 1, time);
    outfile.putVarPart('h',  [0,0,outStep], [Np,Ne,1], h);
    outfile.putVarPart('qx', [0,0,outStep], [Np,Ne,1], qx);
    outfile.putVarPart('qy', [0,0,outStep], [Np,Ne,1], qy);
    outStep = outStep + 1;
    
    contour = contour + 1;
    fprintf('Processing:%f, dt:%f, s:%f...\n', time/FinalTime, dt, s);

end
tElapsed = toc;
outfile.putVarPart('elapsedTime', 0, 1, tElapsed);
%% Assignment
phys.h  = h;
phys.qx = qx;
phys.qy = qy;
% close output file
outfile.CloseFile;
end

function [rk4a, rk4b, rk4c] = RK45coef
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

function s = SWE_PredictWaveSpeed2d(phys, h, Qx, Qy)
% parameters
g        = phys.gra;
minDepth = phys.minht;
dryflag  = h<=minDepth;
dryEle   = any(dryflag);

% estimate of the wave speed
u   = Qx./h; u(:, dryEle) = 0.0;
v   = Qy./h; v(:, dryEle) = 0.0;
spe = sqrt(u.^2 + v.^2);
s   = max(max( spe + sqrt(g*h) ));
end% func

% function F = Fliter(shape, Nc, frac)
% 
% filterdiag = ones(shape.nNode, 1);
% 
% % build exponential filter
% sk = 1; N = shape.nOrder;
% for i=0:N
%   for j=0:N-i
%     if (i+j>=Nc)
%       filterdiag(sk) = frac;
%     end
%     sk = sk+1;
%   end
% end
% 
% F = ( shape.VandMatrix*diag(filterdiag) )/shape.VandMatrix;
% end% func