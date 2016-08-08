function [Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWE_Flux2d(phys, h, qx, qy, dryNodeFlag)
% Function to obtain flux terms for 2 dimensional SWE.
% Input :
%   phys - strucure variable
%   h    - water depth
%   qx   - water flux along x
%   qy   - water flux along y
%   dryNodeFlag - dry flag for each nodes
% 
% Output :
%   F - vector of flux term in x coordinate
%   G - vector of flux term in y coordinate
% 
%% Usages
%   
%   shape        = StdRegions.Triangle(1);
%   [VX,VY,EToV] = Utilities.Mesh.MeshGenTriangle2D(1,-1,1,0);
%   mesh         = MultiRegions.RegionTri(shape, EToV, VX, VY);
%   h  = [0, 1e-4, 1e-3; 1, 1e-4, 1]';
%   qx = [0, 1e-3, 1e-3; 0, 1e-3, 1e-4]';
%   qy = [0, 1e-3, 2e-3; 1, 1e-3, 1e-3]';
%   dryEleFlag = [1, 1, 1; 0, 0, 0]';
% 
%   phys.gra = 9.81;
%   phys.mesh = mesh;
%   [Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWEFlux2d(phys, h, qx, qy, dryEleFlag)
% 

%% Parameters and wet element flag
g           = phys.gra;
minDepth    = phys.minDepth;
wetNodeFlag = h > minDepth;

%% Get the flux terms
% The flux term of SWE is 
% $\mathbf{F}(U) = \left( \begin{array}{c} hu \cr hu^2 + \frac{1}{2} gh^2
% \cr huv \end{array} \right), \quad \mathbf{G}(U) = \left( \begin{array}{c}
%  hv \cr huv \cr hv^2 +\frac{1}{2} gh^2 \end{array} \right)$.
% For dry elements the flux term is zeros.

[r,c] = size(h);
u   = zeros(r,c); v   = zeros(r,c);

% % check any zero depth in wet cells and give a warning
% troubleCell = find(any( h(isWet) <= 0.0 ));
% if troubleCell
%     warning('There exists zero depth in wet cell: %d', troubleCell);
% end% if

% flow rate
u(wetNodeFlag)  = qx(wetNodeFlag)./h(wetNodeFlag);
v(wetNodeFlag)  = qy(wetNodeFlag)./h(wetNodeFlag);

% flux for h
Fh = qx;
Gh = qy;

% flux for hu
Fqx = u.^2.*h + 0.5.*g.*h.^2; 
Gqx = h.*u.*v;

% flux for hv
Fqy = h.*u.*v;
Gqy = v.^2.*h + 0.5.*g.*h.^2; 

% eliminate dry cell flux
Fqx(dryNodeFlag) = 0.0; Gqx(dryNodeFlag) = 0.0;
Fqy(dryNodeFlag) = 0.0; Gqy(dryNodeFlag) = 0.0;

end