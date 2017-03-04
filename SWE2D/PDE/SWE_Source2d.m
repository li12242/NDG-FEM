function [Sh, Sqx, Sqy] = SWE_Source2d(phys, mesh, h, qx, qy, bot)
% Calculation of source terms
% Input: 
%   phys - structure variable, it contains
%       |
%       |   - gra: gravity acceleration
% 
%   mesh - mesh object
%   h    - water depth
%   bot  - bottom elvation
%   dryEleFlag - flag for dry elements
% Output:
%   [Sh, Sqx, Sqy] - source terms for all eqs
% 
%% Usages:
% 
%   shape        = StdRegions.Triangle(1);
%   [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(1,-1,1,0);
%   mesh         = MultiRegions.RegionTri(shape, EToV, VX, VY);
%   bot = mesh.x;
%   h   = [0, 1e-4, 1e-3; 1, 1e-4, 1]';
%   qx  = [0, 1e-3, 1e-3; 0, 1e-3, 1e-4]';
%   qy  = [0, 1e-3, 2e-3; 1, 1e-3, 1e-3]';
%   dryEleFlag = [true, false];
% 
%   phys.gra = 9.81;
%   phys.mesh = mesh;
%   [Sh, Sqx, Sqy] = SWESource2d(phys, mesh, h, bot, dryEleFlag)
% 

% get parameters
g         = phys.gra;
MannCoeff = phys.ManningCoeff;
dryNode   = h<phys.minDepth;

% calculate the gradient of bottom level
Sh  = zeros(size(h));
[Sqx, Sqy] = Grad2D(mesh, bot);

% multiply the depth h and g
Sqx = -g*(h.*Sqx);
Sqy = -g*(h.*Sqy);

qn  = sqrt(qx.^2+qy.^2);
p   = 7/3;
Sqx = Sqx + -g.*MannCoeff.^2.*qx.*qn./(h.^p);
Sqy = Sqy + -g.*MannCoeff.^2.*qy.*qn./(h.^p);

% eleminate dry element terms
Sqx(dryNode) = 0;
Sqy(dryNode) = 0;

end% func