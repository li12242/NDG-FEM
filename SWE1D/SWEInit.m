function physics = SWEInit(physics, nOrder, nEle)
% initialization of the simulation
% Input:
% 
% Output:

caseName = physics.getVal('caseName');
switch caseName
    case 'DamBreakDry'
        FinalTime = 20; % Dam break
        x1 = 0; x2 = 1000; % Dam break
    case 'DamBreakWet'
        FinalTime = 20; % Dam break
        x1 = 0; x2 = 1000; % Dam break
    case 'FlowDump'
        FinalTime = 200; % Flow over dump
        x1 = 0; x2 = 25; % Flow over dump
    case 'ParabolicBowl'
        a = 3000; h0 = 10; g = 9.81;
        T = 2*pi*a/sqrt(2*g*h0);
        FinalTime = 2*T; % Parabolic Bowl
        x1 = -5000; x2 = 5000; % Parabolic Bowl
    case 'LakeAtRest'
        FinalTime = 0.5;
        x1 = 0; x2 = 1;
    case 'TsunamiRunup'
%         FinalTime = 240;
        FinalTime = 10;
        x1 = -500; x2 = 50000;
end% switch
physics.incert('FinalTime', FinalTime);

% set mesh vertex
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nEle);
BC = [2,1; 3,Nv];

if strcmp( caseName, 'TsunamiRunup' )
    % nonuniform mesh discretization
    K1 = floor(K*3/4);
    K2 = floor((K - K1)./2);
    K3 = K - K1 - K2;
    
    xt1 = 2e4; xt2 = 3e4;
    VX(1:(K1+1)) = linspace(x1, xt1, K1+1);
    VX((K1+1):(K1+K2+1)) = linspace(xt1, xt2, K2+1);
    VX((K1+K2+1):(K+1)) = linspace(xt2, x2, K3+1);
end% if

% Initialize solver and construct grid and metric
line = StdRegions.Line(nOrder);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);
physics.incert('mesh', mesh);
physics.incert('VX', VX);
physics.incert('EToV', EToV);

% set bottom topography
bedElva = SetBed(physics);
physics.incert('bedElva', bedElva);

% set initial condition
switch caseName
    case 'DamBreakDry'
        [h, q] = DamBreakInit(mesh, 2);
    case 'DamBreakWet'
        [h, q] = DamBreakInit(mesh, 1);
    case 'FlowDump'
        initCase = 3;
        [h, q] = FlowDumpInit(mesh, bedElva, initCase);
    case 'ParabolicBowl'
        [h, q] = ParaBowlInit(mesh, bedElva);
        
    case 'LakeAtRest'
        [h, q] = LakeAtRestInit(mesh, bedElva);
        
    case 'TsunamiRunup'
        [h, q] = TsunamiRunup(mesh, bedElva);
end% switch
physics.incert('height', h);
physics.incert('flux', q);
end% func

function [h, q] = TsunamiRunup(mesh, bedElva)
q = zeros(size(mesh.x));
data = load('TsunamiRunupInitialCondition.mat');
Interp = griddedInterpolant(data.x, data.eta+5000, 'nearest');
z = Interp(mesh.x(:)); z = reshape(z, size(mesh.x));

h = z - bedElva;
h(h<0) = 0;
end% func

function [h, q] = LakeAtRestInit(mesh, bedElevation)
eta = ones(size(mesh.x));
h = eta - bedElevation;
q = zeros(size(mesh.x));
% correct transition element
transIndex = find(TransiteCellIdentify(mesh, h));
h(h<0) = 0;

np = 10;
[r, w] = Polylib.zwglj(np);
V = StdRegions.Line.getVandMatrix(mesh.Shape.nOrder, r);
for i = 1:numel(transIndex)
    b1 = bedElevation(1, transIndex(i)); b2 = bedElevation(2, transIndex(i));
%     x = (1-r)./2*x1 + (1+r)./2*x2;
    b = (1-r)./2*b1 + (1+r)./2*b2;
    htemp = 1 - b; htemp(htemp<0) = 0;
    temp = htemp.*w;
    hm = (sum( V.* repmat(temp, 1, mesh.Shape.nNode)));
    hm(2:end) = 0;
    h(:, transIndex(i)) = mesh.Shape.VandMatrix*hm';
end% for

end% func

function [h, q] = ParaBowlInit(mesh, bedElva)
q = zeros(size(mesh.x)); %hDelta = 0.0;

a = 3000; h0 = 10; g = 9.81; B = 5; w = sqrt(2*g*h0)./a;
% z = zeros(size(mesh.x));
z = (-2*B.^2 -(4*B*w).*mesh.x)./(4*g);
h = z - bedElva;

% correct transition element
% transIndex = find(TransiteCellIdentify(mesh, h));
h(h<0) = 0;

% np = 10;
% [r, w] = Polylib.zwglj(np);
% V = StdRegions.Line.getVandMatrix(mesh.Shape.nOrder, r);
% for i = 1:numel(transIndex)
%     b1 = bedElva(1, transIndex(i)); b2 = bedElva(end, transIndex(i));
%     z1 = z(1, transIndex(i)); z2 = z(end, transIndex(i));
%     eta = (1-r)./2*z1 + (1+r)./2*z2;
%     b = (1-r)./2*b1 + (1+r)./2*b2;
%     htemp = eta - b; htemp(htemp<0) = 0;
%     temp = htemp.*w;
%     hm = (sum( V.* repmat(temp, 1, mesh.Shape.nNode)));
%     hm(2:end) = 0;
%     h(:, transIndex(i)) = mesh.Shape.VandMatrix*hm';
% end% for

end% func

function [h, q] = FlowDumpInit(mesh, bedElva, initCase)
switch initCase
    case 1 % subcritical flow
        h = 0.5.*ones(size(mesh.x))- bedElva; 
        q = 0.18.*ones(size(mesh.x));
    case 2 % supercritical flow
        h = 2.0.*ones(size(mesh.x))- bedElva; 
        q = 25.0567.*ones(size(mesh.x));
    case 3 % transcritical flow
        h = 0.33.*ones(size(mesh.x))- bedElva;
        q = 0.18.*ones(size(mesh.x));
end% switch
end% func

function [h, q] = DamBreakInit(mesh, initCase)
% Idealized dam break problem of 1D shallow water equation
h = 10.*ones(size(mesh.x)); q = zeros(size(mesh.x));
damPosition = 500;
switch initCase
    case 1 % wet bed
        flag = mesh.x > damPosition;
        h(flag) = 2;
    case 2 % dry bed
        flag = mesh.x > damPosition;
        h(flag) = 1e-4;
end% switch
end% func

function bedElevation = SetBed(physics)
mesh = physics.getVal('mesh');
VX = physics.getVal('VX');
EToV = physics.getVal('EToV');
% Initialize the bed elevation according to the case name
switch physics.getVal('caseName')
    case 'DamBreakDry'
        VB = zeros(size(VX));
%         bedElevation = zeros(size(mesh.x));
    case 'DamBreakWet'
        VB = zeros(size(VX));
%         bedElevation = zeros(size(mesh.x));
    case 'FlowDump'
        VB = zeros(size(VX));
        flag = (VX >= 8) & (VX <=12);
        VB(flag) = 0.2 - 0.05*(VX(flag) -10).^2;        
%         bedElevation = zeros(size(mesh.x));
%         flag = (mesh.x >= 8) & (mesh.x <=12);
%         bedElevation(flag) = 0.2 - 0.05*(mesh.x(flag) -10).^2;
    case 'ParabolicBowl'
        a = 3000; h0 = 10;
        VB = h0.*(VX.^2./a^2 - 1);
%         bedElevation = h0.*(mesh.x.^2./a^2 - 1);
    case 'LakeAtRest'
        VB = zeros(size(VX));
        a = 1.2; rm = 0.4; R = abs(VX - 0.5);
        index = (R < rm);
        VB(index) = a*exp(-0.5./(rm.^2 - R(index).^2))./exp(-0.5./rm^2);
    case 'TsunamiRunup'
        VB = 5000 - 0.1*VX;
end% switch
% all of the bottom level is linear polynomial
physics.incert('VB', VB);
vb = VB(EToV');
bedElevation = 0.5*((1-mesh.Shape.r)*vb(1,:) + (mesh.Shape.r+1)*vb(2,:));
end% func