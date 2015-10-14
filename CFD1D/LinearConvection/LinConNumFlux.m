function cs = LinConNumFlux(mesh, u, c)
% 1D Linear Convection Test Case
% Purpose  : Evaluate numerical flux
% REFERENCE:
% [1] Toro E F. Riemann solvers and numerical methods for fluid dynamics: 
%     a practical introduction[M]. Springer Science & Business Media, 
%     2009. 214-214.

% [^1] L-F
cs = mesh.nx.*u.*(c(mesh.vmapM)+c(mesh.vmapP))./2 ...
    - abs(u*mesh.nx).*(c(mesh.vmapP) - c(mesh.vmapM))./2;
end% func