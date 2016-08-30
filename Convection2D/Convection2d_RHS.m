function rhsVar = Convection2d_RHS(mesh, var, u, v)
% 2D convection problem
% Righ Hand sides

shape = mesh.Shape;

% volume flux
[F, G] = Convection2d_Flux(var, u, v);

% numerical flux
Fs = Convection2d_LF(mesh, var, u, v);
FM = F(mesh.vmapM); 
GM = G(mesh.vmapM);
dF = normal(mesh, FM, GM) - Fs;

% strong form
rhsVar = -( mesh.rx.*(shape.Dr*F) + mesh.sx.*(shape.Ds*F) ...
    + mesh.ry.*(shape.Dr*G) + mesh.sy.*(shape.Ds*G) ) ...
    + shape.LIFT*(dF.*mesh.fScale);

% weak form
% rhsVar = ( mesh.rx.*(tri.Drw*F) + mesh.sx.*(tri.Dsw*F) ...
%     + mesh.ry.*(tri.Drw*G) + mesh.sy.*(tri.Dsw*G) ) ...
%     - (tri.invM*tri.Mes)*(Fs.*mesh.fScale);
end% func

function Fs = Convection2d_LF(mesh, var, u, v)
% L-F flux
uM   = u(mesh.vmapM); 
vM   = v(mesh.vmapM);
varM = var(mesh.vmapM); 
varP = var(mesh.vmapP);
un   = normal(mesh, uM, vM);

% the boundary nodes
bclist = (mesh.vmapM == mesh.vmapP);
varP(bclist) = 0;

[FM, GM] = Convection2d_Flux(varM, uM, vM);
[FP, GP] = Convection2d_Flux(varP, uM, vM);

FsM = normal(mesh, FM, GM);
FsP = normal(mesh, FP, GP);

Fs  = 0.5*(FsM + FsP) - 0.5.*abs(un).*(varP - varM);
end% func

function un = normal(mesh, u, v)
un = u.*mesh.nx + v.*mesh.ny;
end

function [F,G] = Convection2d_Flux(var, u, v)
% Cal flux
% F = u*c; G = v*c
F = var.*u; G = var.*v;
end% func