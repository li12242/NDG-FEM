function phys = SWE_Init1d(phys)
% initialization of the simulation
% Input:
%   phys
% Output:
%   phys

%% Parameters
caseName = phys.name;
N        = phys.n;
Ne       = phys.ne;
phys.spl = [];
%% Init test case
switch caseName
    case 'DamBreakDry'
        [mesh, h, q, bot, ftime, dt, dx] = DamBreakDryInit(phys, N, Ne);
    case 'DamBreakWet'
        [mesh, h, q, bot, ftime, dt, dx] = DamBreakWetInit(N, Ne);
    case 'FlowDump'
        testcase = 1;
        [mesh, h, q, bot, ftime, dt, dx] = FlowDumpInit(N, Ne, testcase);
    case 'ParabolicBowl'
        [mesh, h, q, bot, ftime, dt, dx] = ParabolicBowlInit(phys, N, Ne);
    case 'LakeAtRest'
        [mesh, h, q, bot, ftime, dt, dx] = LakeAtResrInit(N, Ne);
    case 'TsunamiRunup'
        [mesh, h, q, bot, ftime, dt, dx] = TsunamiRunupInit(N, Ne);
    case 'WiderGaussianMound'
        [mesh, h, q, bot, ftime, dt, dx] = WiderGaussianMound(N, Ne);
    case 'LocalGaussianMound'
        [mesh, h, q, bot, ftime, dt, dx] = LocalGaussianMound(N, Ne);
    case 'WiderMovingMound'
        [mesh, h, q, bot, ftime, dt, dx] = WiderMovingMound(N, Ne);
    case 'LocalMovingMound'
        [mesh, h, q, bot, ftime, dt, dx, spl1, spl2] = LocalMovingMound(N, Ne);
        phys.spl = [spl1, spl2]; % add sponge layer object to struct variable
end% switch

%% Assignment
phys.mesh  = mesh;
phys.h     = h;
phys.q     = q;
phys.bot   = bot;
phys.ftime = ftime;
phys.dt    = dt;
phys.dx    = dx;
end% func

function [mesh, h, q, bot, ftime, dt, dx, spl1, spl2] = LocalMovingMound(N, Ne)
% Init mesh
x1                = -640e3; 
x2                =  640e3;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line              = StdRegions.Line(N);
mesh              = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 100.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
R       = 50e3;
eta0    = 1;
x0      = -1358e3;
h       = h + eta0*exp(-(mesh.x-x0).^2/R^2);

% Parameters
ftime   = 29*3600; % 9h
dx      = (x2 - x1)./Ne/N;
dt      = 1;

% sponge layer
% sponge layer #1
xc = mean(mesh.x);
spfilename = 'SWE1D_WiderMovingMound_500_BC1.nc';
xb =  510e3; xe =  640e3;
spBCflag = xc > xb; % sponge layer flags
spl1 = MultiRegions.BoundCondition.SpongeBC1d...
    (mesh, spBCflag, spfilename, xb, xe);
% sponge layer #2
spfilename = 'SWE1D_WiderMovingMound_500_BC2.nc';
xb = -510e3; xe = -640e3;
spBCflag = xc < xb; % sponge layer flags
spl2 = MultiRegions.BoundCondition.SpongeBC1d...
    (mesh, spBCflag, spfilename, xb, xe);
end% func

function [mesh, h, q, bot, ftime, dt, dx] = WiderMovingMound(N, Ne)
% Init mesh
x1                = -1920e3; 
x2                = 1920e3;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line              = StdRegions.Line(N);
mesh              = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 100.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
R       = 50e3;
eta0    = 1;
xc      = -1358e3;
h       = h + eta0*exp(-(mesh.x-xc).^2/R^2);

% Parameters
ftime   = 29*3600; % 9h
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

function [mesh, h, q, bot, ftime, dt, dx] = LocalGaussianMound(N, Ne)
% Init mesh
x1                = -640e3; 
x2                =  640e3;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line              = StdRegions.Line(N);
mesh              = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 100.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
R       = 50e3;
eta0    = 1;
h       = h + eta0*exp(-mesh.x.^2/R^2);

% Parameters
ftime   = 9*3600; % 9h
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

function [mesh, h, q, bot, ftime, dt, dx] = WiderGaussianMound(N, Ne)
% Init mesh
x1                = -1920e3; 
x2                = 1920e3;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line              = StdRegions.Line(N);
mesh              = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 100.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
R       = 50e3;
eta0    = 1;
h       = h + eta0*exp(-mesh.x.^2/R^2);

% Parameters
ftime   = 9*3600; % 9h
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for Tsunami-Runup
function [mesh, h, q, bot, ftime, dt, dx] = TsunamiRunupInit(N, Ne)
% Init mesh
x1                = -500; 
x2                = 50000;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
rm                = 1.5; 
VX                = MeshMapping(rm, Nv, x1, x2);
line              = StdRegions.Line(N);
mesh              = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                 = mesh.Shape.r;
VB                = 5000 - 0.1*VX;
vb                = VB(EToV');
bot               = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
q                 = zeros(size(mesh.x));
data              = load('TsunamiRunupInitialCondition.mat');
Interp            = griddedInterpolant(data.x, data.eta+5000, 'nearest');
z                 = Interp(mesh.x(:)); 
z                 = reshape(z, size(mesh.x));
h                 = z - bot;
h(h<0)            = 0;

% Parameters
ftime   = 360;
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for lake at rest
function [mesh, h, q, bot, ftime, dt, dx] = LakeAtResrInit(N, Ne)

% Init mesh
x1               = 0; 
x2               = 1;
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line             = StdRegions.Line(N);
mesh             = MultiRegions.RegionLine(line, EToV, VX);

% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
a                = 1.2; 
rm               = 0.4; 
R                = abs(VX - 0.5);
index            = (R < rm);
VB(index)        = a*exp(-0.5./(rm.^2 - R(index).^2))./exp(-0.5./rm^2);
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );

% Initial condition
eta              = ones(size(mesh.x));
h                = eta - bot;
q                = zeros(size(mesh.x));
% correct transition element
isWet            = WetDryJudge(mesh, h, physics);
transIndex       = find(TransiteCellIdentify(mesh, isWet));
h(h<0)           = 0;
% reconstruct transition cell
np               = 10; % integral points
[r, w]           = Polylib.zwglj(np);
for i = 1:numel(transIndex)
    ind   = transIndex(i);
    b1    = bot(1,   ind); 
    b2    = bot(end, ind);
    b     = (1-r)./2*b1 + (1+r)./2*b2;
    htemp = 1 - b; htemp(htemp<0) = 0;
    hmean = sum(htemp.*w)/2;
    h(:, ind) = hmean;
end% for

% periodic BC
mesh.vmapP(1)    = mesh.nNode;
mesh.vmapP(end)  = 1;

% Parameters
ftime   = 0.5;
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for flow over dump problem
function [mesh, h, q, bot, ftime, dt, dx] = FlowDumpInit(N, Ne, initCase)
% Input
%   N         - polynomial order
%   Ne        - number of elements
%   InitCase  - different test case
% 

% Init mesh
x1               = 0; 
x2               = 25; % Flow over dump
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line             = StdRegions.Line(N);
mesh             = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
flag             = (VX >= 8)&(VX <=12);
VB(flag)         = 0.2 - 0.05*(VX(flag) -10).^2;  
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Init condition
switch initCase
    case 1 % subcritical flow
        h = 0.5.*ones(size(mesh.x))- bot; 
        q = 0.18.*ones(size(mesh.x));
    case 2 % supercritical flow
        h = 2.0.*ones(size(mesh.x))- bot; 
        q = 25.0567.*ones(size(mesh.x));
    case 3 % transcritical flow
        h = 0.33.*ones(size(mesh.x))- bot;
        q = 0.18.*ones(size(mesh.x));
end% switch
% Parameters
ftime   = 200; % Flow over dump
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for parabolic bowl oscillation flow problem
function [mesh, h, q, bot, ftime, dt, dx] = ParabolicBowlInit(phys, N, Ne)
% Init mesh for parabolic bowl problem
a                = 3000; 
h0               = 10; 
g                = phys.gra;
T                = 2*pi*a/sqrt(2*g*h0);
x1               = -5000; 
x2               = 5000; % Parabolic Bowl
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line             = StdRegions.Line(N);
mesh             = MultiRegions.RegionLine(line, EToV, VX);
% Init bottom level
r                = mesh.Shape.r;
VB               = h0.*(VX.^2./a^2 - 1);
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Init Condition
q = zeros(size(mesh.x)); %hDelta = 0.0;
B = 5; 
w = sqrt(2*g*h0)./a;
z = (-2*B.^2 -(4*B*w).*mesh.x)./(4*g);
h = z - bot;
% correct transition element
h(h<0) = 0;

% Parameters
ftime = 2*T; % Parabolic Bowl
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for dry dam break problem
function [mesh, h, q, bot, ftime, dt, dx] = DamBreakDryInit(phys, N, Ne)
% Init mesh for Dam break
x1               = 0; 
x2               = 1000; 
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line             = StdRegions.Line(N);
mesh             = MultiRegions.RegionLine(line, EToV, VX);

% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 10.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
dampos  = 500;
xc      = mean(mesh.x);
h(:, xc > dampos) = phys.minDepth;

% Parameters
ftime   = 20;
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% Init test case for wet dam break problem
function [mesh, h, q, bot, ftime, dt, dx] = DamBreakWetInit(N, Ne)
% Init mesh for Dam break
x1               = 0; 
x2               = 1000; 
[~, VX, ~, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, Ne);
line             = StdRegions.Line(N);
mesh             = MultiRegions.RegionLine(line, EToV, VX);

% Init bottom level
r                = mesh.Shape.r;
VB               = zeros(size(VX));
vb               = VB(EToV');
bot              = 0.5*( (1-r)*vb(1,:) + (r+1)*vb(2,:) );
% Initial condition
h       = 10.*ones(size(mesh.x)); 
q       = zeros(size(mesh.x));
dampos  = 500;
xc      = mean(mesh.x);
h(:, xc > dampos) = 2;

% Parameters
ftime   = 20;
dx      = (x2 - x1)./Ne/N;
dt      = 1;
end% func

%% MeshMapping
function x = MeshMapping(rm, Nv, x1, x2)
fp = @(r) exp(r) - 1;
fm = @(r) 1 - exp(-r);
xm = @(y) - log(1-y);
% fp = @(r) r./3;
% fm = @(r) r./3;
% xm = @(y) 3.*y;

rmin = xm(fp(rm)/x2*x1);
r = linspace(rmin, rm, Nv);
[~, index] = min(abs(r)); r(index) = 0;
x = zeros(size(r));
x(r > 0) = fp( r(r>0) )*x2/fp(r(end));
x(r < 0) = fm( r(r<0) )*x1/fm(r(1));
end% func