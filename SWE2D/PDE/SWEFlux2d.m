function [Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWEFlux2d(phys, h, qx, qy)
% Function to obtain flux terms for 2 dimensional SWE.
% Input :
%   phys - strucure variable
%   h    - water depth
%   qx   - water flux along x
%   qy   - water flux along y
% 
% Output :
%   F - vector of flux term in x coordinate
%   G - vector of flux term in y coordinate
% 

%% Parameters and dry element flag
% 
g = phys.gra;
dryEleFlag = IsDry(mesh, h, minDepth);
isWet = ~dryEleFlag;

%% Get the flux terms
% The flux term of SWE is 
% $\mathbf{F}(U) = \left( \begin{array}{c} hu \cr hu^2 + \frac{1}{2} gh^2
% \cr huv \end{array} \right), \quad \mathbf{G}(U) = \left( \begin{array}{c}
%  hv \cr huv \cr hv^2 +\frac{1}{2} gh^2 \end{array} \right)$.
% For dry elements the flux term is zeros.
% 

[r,c] = size(h);
Fqx = zeros(r,c); Fqy = zeros(r,c);
Gqx = zeros(r,c); Gqy = zeros(r,c);

u  = qx./h; v = qy./h;
u(:, dryEleFlag) = 0;
v(:, dryEleFlag) = 0;

% flux for h
Fh = h.*u;
Gh = h.*v;

% flux for hu
Fqx(:, isWet) = qx(:, isWet).^2./h(:, isWet) + 0.5.*g.*h(:, isWet).^2; 

Gqx(:, isWet) = qx(:, isWet).*qy(:, isWet)./h(:, isWet);

% flux for hv
Fqy(:, isWet) = qx(:, isWet).*qy(:, isWet)./h(:, isWet);
Gqy(:, isWet) = qy(:, isWet).^2./h(:, isWet) + 0.5.*g.*h(:, isWet).^2; 

end