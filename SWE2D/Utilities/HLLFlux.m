function [Fhs, Fqxs, Fqys] = HLLFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryEleFlag)
% Calculate the HLL flux function
% Input
%   hM    - water depth in local element 
%   hP    - water depth in adjacent element
%   QxM   - water flux along x coordinate
%   QxP   - water flux along y coordinate
%   QyM   - 
%   QyP   -
%   dryM  - dry node flag for local nodes
%   dryP  - dry node flag for adjacent nodes
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

gra      = phys.gra;
minDepth = phys.minDepth;
dryM     = (hM <= minDepth);
dryP     = (hP <= minDepth);


%% Estimate wave speed
% The wave speed is estimated from Fraccarollo and Toro (1995) basing on
% the Equation
% 
% $S_L = min(u^- - \sqrt{gh^-}, u^* - c^*) \quad 
% S_R = max(u^+ + \sqrt{gh^+}, u^* + c^*)$
% 
% For the dry bed on the right side of a node ($h^- > 0$ and $h^+=0$) the
% wave speed are given by
% 
% $S_L = u^- - \sqrt{gh^-} \quad S_R = u^- + \sqrt{2gh^-}$
% 
% For the dry bed on the left side of a node ($h^- > 0$ and $h^+=0$) the
% wave speed are given by
% 
% $S_L = u^+ - \sqrt{2gh^+} \quad S_R = u^+ + \sqrt{gh^+}$
% 

uM = QxM./hM; uM(dryM) = 0.0;
uP = QxP./hP; uP(dryP) = 0.0;

us = 0.5*( uM + uP ) + ( sqrt(gra*hM) - sqrt(gra*hP) );
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
SP(flag) = uP(flag) +   sqrt(gra.*hP(flag) );

% right is dry
flag     = ~dryM & dryP; 
SM(flag) = uM(flag) -   sqrt(gra.*hM(flag) );
SP(flag) = uM(flag) + 2*sqrt(gra.*hM(flag) );

%% Calculate the numerical flux
% The HLL flux function is given by 
% 
% $F^{HLL} = \left\{ \begin{array}{cc}
% F^- & S_L \geq 0 \cr 
% \frac{S_R F^{-} - S_L F^{+} + S_L S_R(U^+ - U^-)}{S_R S_L} & S_L < 0 < S_R \cr
% F^+  & S_R \leq 0 \end{array} \right.$
% 
% 
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


flag      = (SM < 0) & (SP > 0);
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
flag       = dryM & dryP;
Fhs(flag)  = 0.0;
Fqxs(flag) = 0.0;
Fqys(flag) = 0.0;

end% func