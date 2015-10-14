function rhsVar = Convection2DRHS(mesh, var, time, Speed)
% 2D convection problem
% Righ Hand sides

tri = mesh.Shape;

Fs = ConvectionLF(mesh, var, Speed);
Fs = BoundaryCondition(mesh, Fs, time, Speed);

[F, G] = ConvectionFlux(var, Speed);

% strong form

% [FM, GM] = ConvectionFlux(var(mesh.vmapM), Speed);
% dF = normal(mesh, FM, GM) - Fs;
% dF = lf(mesh, var, Speed);

% rhsVar = -( mesh.rx.*(tri.Dr*F) + mesh.sx.*(tri.Ds*F) ...
%     + mesh.ry.*(tri.Dr*G) + mesh.sy.*(tri.Ds*G) ) ...
%     + tri.invM*tri.Mes*(dF.*mesh.fScale);

% weak form
rhsVar = ( mesh.rx.*(tri.Drw*F) + mesh.sx.*(tri.Dsw*F) ...
    + mesh.ry.*(tri.Drw*G) + mesh.sy.*(tri.Dsw*G) ) ...
    - (tri.invM*tri.Mes)*(Fs.*mesh.fScale);
end% func

function Fs = ConvectionLF(mesh, var, Speed)
% L-F flux
un = normal(mesh, Speed(1), Speed(2));
[FM, GM] = ConvectionFlux(var(mesh.vmapM), Speed);
[FP, GP] = ConvectionFlux(var(mesh.vmapP), Speed);
FsM = normal(mesh, FM, GM); 
FsP = normal(mesh, FP, GP); 
varM = var(mesh.vmapM); varP = var(mesh.vmapP);
Fs = 0.5*(FsM + FsP) - 0.5.*abs(un).*(varP - varM);
end% func

function Fs = ConvectionCentral(mesh, var, Speed)
% central flux
[FM, GM] = ConvectionFlux(var(mesh.vmapM), Speed);
[FP, GP] = ConvectionFlux(var(mesh.vmapP), Speed);
FsM = normal(mesh, FM, GM); 
FsP = normal(mesh, FP, GP); 

Fs = 0.5*(FsM + FsP);
end% func

function Fs = BoundaryCondition(mesh, Fs, time, Speed)
% Inflow BC
% c = sin(t);
c = sin(- pi*time);%*sin(2*pi*mesh.y(mesh.vmapM(mesh.mapI)));
[F,G] = ConvectionFlux(c, Speed);
Fs(mesh.mapI) = F.*mesh.nx(mesh.mapI) + G.*mesh.ny(mesh.mapI);
end% func


function un = normal(mesh, u, v)
un = u.*mesh.nx + v.*mesh.ny;
end

function [F,G] = ConvectionFlux(var, Speed)
% Cal flux
% F = u*c; G = v*c
F = var*Speed(1); G = var*Speed(2);
end% func