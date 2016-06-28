function [Fhs, Fqxs, Fqys] = HLLFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryM, dryP)
% Calculate the HLL flux function
% Input:
%   hM - water depth in local element 
%   hP - water depth in adjacent element
% 
% 
% Output:
% 
%% Usages
% 
%   shape        = StdRegions.Triangle(1);
%   [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(1,-1,1,0);
%   mesh         = MultiRegions.RegionTri(shape, EToV, VX, VY);
%   h   = [0, 1e-4, 1e-3; 1, 1e-4, 1]';
%   qx  = [0, 1e-3, 1e-3; 0, 1e-3, 1e-4]';
%   qy  = [0, 1e-3, 2e-3; 1, 1e-3, 1e-3]';
%   dryEleFlag  = [true, false];
%   dryNodeFlag = repmat(dryEleFlag, mesh.Shape.nNode, 1);
%   dryM        = dryNodeFlag(mesh.vmapM);
%   dryP        = dryNodeFlag(mesh.vmapP);
% 
%   QxM =  qx(mesh.vmapM).*mesh.nx + qy(mesh.vmapM).*mesh.ny;
%   QyM = -qx(mesh.vmapM).*mesh.ny + qy(mesh.vmapM).*mesh.nx;
% 
%   QxP =  qx(mesh.vmapP).*mesh.nx + qy(mesh.vmapP).*mesh.ny;
%   QyP = -qx(mesh.vmapP).*mesh.ny + qy(mesh.vmapP).*mesh.nx;
% 
%   hM  = h(mesh.vmapM);
%   hP  = h(mesh.vmapP);
% 
%   phys.gra = 9.81;
%   phys.mesh = mesh;
%   [Fhs, Fqxs, Fqys] = HLLFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryM, dryP)
% 

gra  = phys.gra;

%% Estimate wave speed
uM = QxM./hM; uM(dryM) = 0.0;
uP = QxP./hP; uP(dryP) = 0.0;

us = 0.5*(uM + uP) + (sqrt(gra*hM) - sqrt(gra*hP));
cs = 0.5*( sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(uM - uP);

% allocate memory
SM = zeros(size(uM)); 
SP = zeros(size(uP));

% both is wet
flag     = ~(dryM | dryP);
SM(flag) = uM(flag) - sqrt(gra.*hM(flag));
SM(flag) = min(SM(flag), us(flag) - cs(flag) );

SP(flag) = uP(flag) + sqrt(gra.*hP(flag));
SP(flag) = max(SP(flag), us(flag) + cs(flag) );

% left is dry
flag     = dryM & ~dryP;
SM(flag) = uP(flag) - 2*sqrt(gra.*hP(flag) ); 
SP(flag) = uP(flag) + sqrt(gra.*hP(flag) );

% right is dry
flag     = ~dryM & dryP;
SP(flag) = uM(flag) + 2*sqrt(gra.*hM(flag) ); 
SM(flag) = uM(flag) - sqrt(gra.*hM(flag) );

%% Calculate the numerical flux
[FhM, FqxM, FqyM, ~, ~, ~] = SWEFlux2d(phys, hM, QxM, QyM, dryM);
[FhP, FqxP, FqyP, ~, ~, ~] = SWEFlux2d(phys, hP, QxP, QyP, dryP);

Fhs  = zeros(size(FhM)); 
Fqxs = zeros(size(FqxM));
Fqys = zeros(size(FqyM));

flag       = (SM >= 0);
Fhs (flag) = FhM (flag);
Fqxs(flag) = FqxM(flag);
Fqys(flag) = FqyM(flag);


flag       = (SP <= 0);
Fhs (flag) = FhP (flag);
Fqxs(flag) = FqxP(flag);
Fqys(flag) = FqyP(flag);


flag = (SM < 0) & (SP > 0);
Fhs(flag) = (- SM(flag).*FhP(flag)+ SP(flag).*FhM(flag)...
    + SP(flag).*SM(flag).*( hP(flag) - hM(flag) ) )...
    ./(SP(flag) - SM(flag));
Fqxs(flag) = (- SM(flag).*FqxP(flag)+ SP(flag).*FqxM(flag)...
    + SP(flag).*SM(flag).*( QxP(flag) - QxM(flag) ) )...
    ./(SP(flag) - SM(flag));
Fqys(flag) = (- SM(flag).*FqyP(flag)+ SP(flag).*FqyM(flag)...
    + SP(flag).*SM(flag).*( QyP(flag) - QyM(flag) ) )...
    ./(SP(flag) - SM(flag));

% Eliminate flux between dry elements
dryInd       = dryM & dryP;
Fhs(dryInd)  = 0.0;
Fqxs(dryInd) = 0.0;
Fqys(dryInd) = 0.0;

end% func