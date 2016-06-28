function [Fhs, Fqxs, Fqys] = LLFFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryM, dryP)
% Input:
%   h  - water depth
%   Qx - water flux
%   Qy - water flux
% 
% Output:
%   Fhs  - 
%   Fqxs -
%   Fqys - 
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
%   [Fhs, Fqxs, Fqys] = LLFFlux(phys, hM, hP, QxM, QxP, QyM, QyP, dryM, dryP)
% 

g  = phys.gra;

%% Estimate wave speed
uM = QxM./hM; uM(dryM) = 0.0;
uP = QxP./hP; uP(dryP) = 0.0;

SM = abs(uM) + sqrt(g*hM); 
SP = abs(uP) + sqrt(g*hP);
% for dry elements, the local wave speed is zero
SM(dryM) = 0; 
SP(dryP) = 0;
S = max(SM, SP);

%% The Lax-Friedrichs flux
% The function is to calculate the flux between two adjacent elements $K^-$
% and $K^+$ given by
% 
% $$\hat{\mathbf{F}} = \frac{F(q_h^-) + F(q_h^+)}{2} + \frac{C}{2}n^{\mp}
% \left( q_h^{\mp} - q_h^{\pm} \right)$$
% 
% where $q_h^-$ and $q_h^+$ are respectively the solutions at common edge of
% element $K^-$ and $K^+$, $n^- = -n^+$ are the outer  normal vector, 
% and $C$ corresponds to the largest value of the absolute maximum eigenvalue
% of the normal flux Jacobian matrix along the common edge, which is given
% by
% 
% $$C = max_{s\in[q_h^-, q_h^+]} \left| \lambda \left( n_x \frac{\partial f}
% {\partial q} +n_y \frac{\partial g}{\partial q} \right) \right| = 
% max_{s\in[q_h^-, q_h^+]} \left( \left| n\cdot u \right| + \left| \sqrt{gH}
% \right| \right)$$
% 
% where $\lambda$ denotes the eigenvalue of the matrix.
% 

[FhM, FqxM, FqyM, ~, ~, ~] = SWEFlux2d(phys, hM, QxM, QyM, dryM);
[FhP, FqxP, FqyP, ~, ~, ~] = SWEFlux2d(phys, hP, QxP, QyP, dryP);

Fhs  = 0.5*(FhM  + FhP  + S.*(hM  - hP ));
Fqxs = 0.5*(FqxM + FqxP + S.*(QxM - QxP));
Fqys = 0.5*(FqyM + FqyP + S.*(QyM - QyP));

% Eliminate flux between dry elements
dryInd       = dryM & dryP;
Fhs(dryInd)  = 0.0;
Fqxs(dryInd) = 0.0;
Fqys(dryInd) = 0.0;

end% func