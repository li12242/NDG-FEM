function [rhsH, rhsQx, rhsQy] = SWE_RHS2d(phys, mesh, h, qx, qy)
% Calculation the RHS of SWE

% parameters
minDepth = phys.minDepth;
shape    = mesh.Shape;
bot      = phys.bot;
% find the dry nodes
dryEleFlag  = SWE_DryEle2d(mesh, h, minDepth);

% flux
[Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWE_Flux2d(phys, h, qx, qy);
% source terms
[Sh, Sqx, Sqy] = SWE_Source2d(phys, mesh, h, bot, dryEleFlag);

% numerical flux
[Fhs, Fqxs, Fqys] = SWE_NumFlux2d(phys, mesh, h, qx, qy);

% Boundary condition
[Fhs, Fqxs, Fqys] = SWE_BC2d(phys, mesh, h, qx, qy, Fhs, Fqxs, Fqys);

% the flux difference
dFh   = Fh(mesh.vmapM).*mesh.nx + Gh(mesh.vmapM).*mesh.ny   - Fhs;
dFqx  = Fqx(mesh.vmapM).*mesh.nx + Gqx(mesh.vmapM).*mesh.ny - Fqxs;
dFqy  = Fqy(mesh.vmapM).*mesh.nx + Gqy(mesh.vmapM).*mesh.ny - Fqys;

%% strong form
divFh  = Div2D(mesh, Fh, Gh);
rhsH   = - divFh + Sh + shape.LIFT*(dFh.*mesh.fScale);

divFh  = Div2D(mesh, Fqx, Gqx);
rhsQx   = - divFh + Sqx + shape.LIFT*(dFqx.*mesh.fScale);

divFh  = Div2D(mesh, Fqy, Gqy);
rhsQy   = - divFh + Sqy + shape.LIFT*(dFqy.*mesh.fScale);

%% weak form
% divFh  = DivWeak2D(mesh, Fh, Gh);
% rhsH   = divFh + Sh - shape.LIFT*(Fhs.*mesh.fScale);
% 
% divFh  = DivWeak2D(mesh, Fqx, Gqx);
% rhsQx  = divFh + Sqx - shape.LIFT*(Fqxs.*mesh.fScale);
% 
% divFh  = DivWeak2D(mesh, Fqy, Gqy);
% rhsQy  = divFh + Sqy - shape.LIFT*(Fqys.*mesh.fScale);
end

function [Fhs, Fqxs, Fqys] = SWE_BC2d(phys, mesh, h, qx, qy, Fhs, Fqxs, Fqys)
%% Parameters
gra  = phys.gra;
hmin = phys.minDepth;

%% Wall
nodeWallM = mesh.vmapM(mesh.mapW);

nx  = mesh.nx(mesh.mapW);
ny  = mesh.ny(mesh.mapW);
% no-slip wall condition
hM  = h(nodeWallM);   hP  = hM;
qxM = qx(nodeWallM);  qxP = -qxM;
qyM = qy(nodeWallM);  qyP = -qyM;

[Fh, Fqx, Fqy] = SWE_Mex_BC2d...
    (hmin, gra, hM, hP, qxM, qxP, qyM, qyP, nx, ny);
% assignment to numerical flux
Fhs(mesh.mapW)  = Fh;
Fqxs(mesh.mapW) = Fqx;
Fqys(mesh.mapW) = Fqy;
%% Outflow

%% Inflow
end% func