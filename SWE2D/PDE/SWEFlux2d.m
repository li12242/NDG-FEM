function [Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWEFlux2d(phys, h, qx, qy, dryNodeFlag)
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
g          = phys.gra;
minDepth   = phys.minDepth;
dryNodeFlag(h <= minDepth) = 1;

%% Get the flux terms
% The flux term of SWE is 
% $\mathbf{F}(U) = \left( \begin{array}{c} hu \cr hu^2 + \frac{1}{2} gh^2
% \cr huv \end{array} \right), \quad \mathbf{G}(U) = \left( \begin{array}{c}
%  hv \cr huv \cr hv^2 +\frac{1}{2} gh^2 \end{array} \right)$.
% For dry elements the flux term is zeros.

[r,c] = size(h);
Fh  = zeros(r,c); Gh  = zeros(r,c);
Fqx = zeros(r,c); Fqy = zeros(r,c);
Gqx = zeros(r,c); Gqy = zeros(r,c);

u  = qx./h; v = qy./h;
u(dryNodeFlag) = 0;
v(dryNodeFlag) = 0;

% check any zero depth in wet cells and give a warning
isWet       = ~dryNodeFlag;
troubleCell = find(any( h(isWet) <= 0.0 ));
if troubleCell
    warning('There exists zero depth in wet cell: %d', troubleCell);
end% if

% flux for h
Fh(isWet) = h(isWet).*u(isWet);
Gh(isWet) = h(isWet).*v(isWet);

% flux for hu
Fqx(isWet) = qx(isWet).^2./h(isWet) + 0.5.*g.*h(isWet).^2; 
Gqx(isWet) = qx(isWet).*qy(isWet)./h(isWet);

% flux for hv
Fqy(isWet) = qx(isWet).*qy(isWet)./h(isWet);
Gqy(isWet) = qy(isWet).^2./h(isWet) + 0.5.*g.*h(isWet).^2; 

end