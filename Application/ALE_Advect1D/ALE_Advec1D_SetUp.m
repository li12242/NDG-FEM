function ALE_Advec1D_SetUp(nOrder, nEle, filename)
% Advction of ALE form Equation
% 
% Reference: 
% [1]: Ren, X., 2015. A multi?dimensional high?order DG?ALE method based on
%       gas?kinetic theory with application to oscillating airfoils 
%       - Google Ñ§ÊõËÑË÷ [WWW Document]. URL http://scholar.glgoo.org/scholar?hl=zh-CN&q=A+multi%E2%80%90dimensional+high%E2%80%90order+DG%E2%80%90ALE+method+based+on+gas%E2%80%90kinetic+theory+with+application+to+oscillating+airfoils&btnG=&lr= (accessed 12.5.15).
% [2]: 

physics = Utilities.varGroup; 
caseName = 'ALE-Advect1D'; physics.incert('caseName', caseName);

%% set domain

FinalTime = 10; physics.incert('FinalTime', FinalTime); % time domain

x1 = 0; x2 = 2*pi; % space domain
N = nOrder; nElement = nEle; delta_x = (x2 - x1)/nElement;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1; 3,Nv];

line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
physics.incert('mesh', mesh);

% initialization
u = sin(mesh.x); a = ones(size(mesh.x));
physics.incert('speed', a);
physics.incert('scalar', u);

%% set output file
fileName = [filename,'.nc'];
ncfile = CreateOutputFile(fileName, mesh);


%% time evolution
[rk4a, rk4b, rk4c] = getRK_coefficient;

% set time interval
min_x = min( abs(mesh.x(2) - mesh.x(1)) ); CFL = 0.5;
dt = CFL*(min_x./a(1));
nt = floor(FinalTime/dt);
time = 0;

for it = 1:nt
    
    u_temp = u.*mesh.J;
    resu = zeros(size(u));
    
    % cal new mesh & velocity
    % mesh update
    mesh_new = mesh_deform(mesh, time+dt, FinalTime, VX, delta_x, EToV);
    % fix mesh
%     mesh_new = mesh;
    if dt > 0
        vg = (mesh_new.x - mesh.x)/dt;           
        at = a - vg;
    else 
        at = a;
    end%
%     mesh_temp = mesh_new;
    for INTRK = 1:4
        timelocal = time + rk4c(INTRK)*dt;
        mesh_temp = mesh_ratio(mesh, mesh_new, rk4c(INTRK));

%         cal RHS
        [rhsu] = ALE_AdvecRHS1D(mesh_temp, u_temp./mesh_temp.J, timelocal, at);
        rhsu = rhsu.*mesh_temp.J;
        resu = resu + dt*rk4b(INTRK)*rhsu;
        u_temp = u.*mesh.J + dt*rk4a(INTRK)*rhsu;
    end;
    
    u = u.*mesh.J./mesh_new.J + resu./mesh_new.J;
    mesh = mesh_new;
    
    subplot(2,1,1)
    plot(mesh.x, u, 'b.', mesh.x, sin(mesh.x - time -dt),'r'); 
    set(gca, 'YLim', [-1.1, 1.1]); drawnow;
    subplot(2,1,2)
    plot(mesh.x, zeros(size(mesh.x)), 'o')
    % Increment time
    time = time+dt;
    % 
    StoreVar(fileName, ncfile, mesh.x, u, time, it-1);
    fprintf('Processing: %f...\n', it/nt)
    
end% for

end% func

function [rhsu] = ALE_AdvecRHS1D(mesh, u, timelocal, a)

line = mesh.Shape; vmapI = 1;

% form field differences at faces
alpha=0;
% du = zeros(size(mesh.nx));
du = (u(mesh.vmapM)-u(mesh.vmapP)).*...
    (a(mesh.vmapM).*mesh.nx-(1-alpha)*abs(a(mesh.vmapM).*mesh.nx))/2;

% impose boundary condition at x=0
uin = -sin(a(1)*timelocal);

% uin = 0;
du (mesh.mapI) = (u(vmapI)- uin ).*(a(1)*mesh.nx(mesh.mapI)...
    -(1-alpha)*abs(a(1)*mesh.nx(mesh.mapI)))/2;
du (mesh.mapO) = 0;

% compute right hand sides of the semi-discrete PDE
rhsu = -mesh.rx.*(line.Dr*(a.*u)) + line.invM*line.Mes*( du.*mesh.fScale );
end% func

function mesh_temp = mesh_ratio(mesh, new_mesh, p)
% update mesh object with linear interpolation:
% 1. coordinate position
% 2. Jacobian scals of surface integral
% 
% Input:
%   mesh        - mesh object
%   new_mesh    - new mesh
%   p           - ratio
% Output:
%   mesh_temp   - interpolate mesh object
mesh_temp = mesh;
vx = mesh.x(mesh.Shape.getFaceListToNodeList, :);
vx_new = new_mesh.x(new_mesh.Shape.getFaceListToNodeList, :);
vx_temp = (1 - p)*vx + p*vx_new;

[mesh_temp.x, mesh_temp.rx, mesh_temp.J] = mesh.Shape.getEleGeometric(vx_temp);
mesh_temp.fScale = mesh_temp.sJ./mesh_temp.J(mesh_temp.vmapM);

end% func

function mesh_out = mesh_deform(mesh, time, FinalTime, VX, delta_x, EToV)
% get new mesh object, update:
% 1. coordinate position
% 2. Jacobian scals of surface integral
% 
% Input:
%   mesh        - mesh object
%   time        - local time
%   FinalTime   - simulation end time, use to calculate period
%   VX          - original position
%   delta_x     - original element size, use to calculate amplitude
% Output:
%   mesh        - new mesh object
% 
mesh_out = mesh;
w = 2*pi/(FinalTime/3); % period T = finalTime/3 = 3.3s
a = delta_x/4; % amplitude

temp = VX + a.*sin( w*time + VX );
VX(2:end-1) = temp(2:end-1); % No 1st vertice doesn't move
vx = zeros(2, mesh.nElement);
vx(1,:) = VX(EToV(:,1)); vx(2,:) = VX(EToV(:,2));

[mesh_out.x, mesh_out.rx, mesh_out.J] = mesh.Shape.getEleGeometric(vx);
mesh_out.fScale = mesh_out.sJ./mesh_out.J(mesh_out.vmapM);

end% func

function [rk4a, rk4b, rk4c] = getRK_coefficient
rk4c = [0, 1/2, 1/2, 1];
rk4a = [1/2, 1/2, 1, 1];
rk4b = [1/6, 1/3, 1/3, 1/6];

%% LSERK
% rk4a = [            0.0 ...
%         -567301805773.0/1357537059087.0 ...
%         -2404267990393.0/2016746695238.0 ...
%         -3550918686646.0/2091501179385.0  ...
%         -1275806237668.0/842570457699.0];
% rk4b = [ 1432997174477.0/9575080441755.0 ...
%          5161836677717.0/13612068292357.0 ...
%          1720146321549.0/2090206949498.0  ...
%          3134564353537.0/4481467310338.0  ...
%          2277821191437.0/14882151754819.0];
% rk4c = [             0.0  ...
%          1432997174477.0/9575080441755.0 ...
%          2526269341429.0/6820363962896.0 ...
%          2006345519317.0/3224310063776.0 ...
%          2802321613138.0/2924317926251.0];
end