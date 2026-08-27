function warp = warpfactor( N, rout )
%> @brief Compute the scaled warp function at order N for the given points.
%> @details Shared by the triangle node placement (see ndgcell.triNodeCoor),
%> after Hesthaven & Warburton (Nodes2D.m).
% ======================================================================
%> This function is part of the NDG-FEM software.
% ======================================================================

% Compute LGL and equidistant node distribution
[LGLr,~] = zwglj(N+1); req = linspace(-1,1,N+1)';

% Compute V based on req
Veq = ndgcell.lineVandMatrix(N, req);
% Evaluate Lagrange polynomial at rout
Nr = length(rout);
Pmat = zeros(N+1,Nr);
for i=1:N+1
  Pmat(i,:) = JacobiP(rout, 0, 0, i-1);
end
Lmat = Veq'\Pmat;
% Compute warp factor
warp = Lmat'*(LGLr - req);
% Scale factor
zerof = (abs(rout)<1.0-1.0e-10); sf = 1.0 - (zerof.*rout).^2;
warp = warp./sf;
end% func
