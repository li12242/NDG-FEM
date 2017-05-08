function [rhsT, p, q] = Convection2d_RHS(mesh, T, u, v, Dx, Dy)
% 2D convection problem
% Righ Hand sides

shape = mesh.Shape;

dT = (T(mesh.vmapM)-T(mesh.vmapP))/2.0;
p = Dx.*( mesh.rx.*(shape.Dr*T) + mesh.sx.*(shape.Ds*T)...
    - shape.LIFT*(mesh.fScale.*(mesh.nx.*dT)) );
q = Dy.*( mesh.ry.*(shape.Dr*T) + mesh.sy.*(shape.Ds*T)...
    - shape.LIFT*(mesh.fScale.*(mesh.ny.*dT)) );

dp = (p(mesh.vmapM)-p(mesh.vmapP))/2.0;
dq = (q(mesh.vmapM)-q(mesh.vmapP))/2.0;
dps = mesh.nx.*dp + mesh.ny.*dq;
% flux term
[F, G] = Convection2d_Flux(T, u, v);

% numerical flux
Fs = Convection2d_LF(mesh, T, u, v);
FM = F(mesh.vmapM); GM = G(mesh.vmapM);
dF = normal(mesh, FM, GM) - Fs;

% strong form
rhsT = -( ...
    + mesh.rx.*(shape.Dr*F) + mesh.sx.*(shape.Ds*F) ...
    + mesh.ry.*(shape.Dr*G) + mesh.sy.*(shape.Ds*G) ) ...
    + shape.LIFT*(dF.*mesh.fScale);
    
rhsT = rhsT + mesh.rx.*(shape.Dr*p) + mesh.sx.*(shape.Ds*p)...
    + mesh.ry.*(shape.Dr*q) + mesh.sy.*(shape.Ds*q)...
    - shape.LIFT*(mesh.fScale.*dps);
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