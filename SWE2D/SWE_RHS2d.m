function [rhsH, rhsQx, rhsQy] = SWE_RHS2d(phys, mesh, h, qx, qy, time)
% Calculation the RHS of SWE

% parameters
minDepth = phys.minDepth;
shape    = mesh.Shape;
bot      = phys.bot;

% find the dry nodes
% dryEleFlag  = SWE_DryEle2d(mesh, h, minDepth);

% flux
[Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWE_Flux2d(phys, h, qx, qy);
% source terms
[Sh, Sqx, Sqy] = SWE_Source2d(phys, mesh, h, qx, qy, bot);

% numerical flux
[Fhs, Fqxs, Fqys] = SWE_NumFlux2d(phys, mesh, h, qx, qy);

% Boundary condition
[Fhs, Fqxs, Fqys] = SWE_BC2d...
    (phys, mesh, h, qx, qy, bot, Fhs, Fqxs, Fqys, time);

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

function [Fhs, Fqxs, Fqys] = SWE_BC2d...
    (phys, mesh, h, qx, qy, bot, Fhs, Fqxs, Fqys, time)
%% Parameters
gra  = phys.gra;
hmin = phys.minDepth;

%% Wall
nodeWallM = mesh.vmapM(mesh.mapW);

nx  = mesh.nx(mesh.mapW);
ny  = mesh.ny(mesh.mapW);

hM  = h(nodeWallM);   hP  = hM;
% no-slip wall condition
% qxM = qx(nodeWallM);  qxP = -qxM;
% qyM = qy(nodeWallM);  qyP = -qyM;

% slip condition
qxM = qx(nodeWallM); qyM = qy(nodeWallM);
qnM =  0; % outward normal flux
qvM = -qxM.*ny + qyM.*nx; % outward tangential flux
qxP = (-qnM).*nx - qvM.*ny;
qyP = (-qnM).*ny + qvM.*nx;

[Fh, Fqx, Fqy] = SWE_Mex_BC2d...
    (hmin, gra, hM, hP, qxM, qxP, qyM, qyP, nx, ny);
% assignment to numerical flux
Fhs(mesh.mapW)  = Fh;
Fqxs(mesh.mapW) = Fqx;
Fqys(mesh.mapW) = Fqy;
%% Outflow

%% Inflow
nodeInM = mesh.vmapM(mesh.mapI);
if (strncmp(phys.casename, 'ObliqueHydraulicJump', 20))
    nx  = mesh.nx(mesh.mapI);
    ny  = mesh.ny(mesh.mapI);
    hM  = h(nodeInM);   hP  = 2.*phys.hin - hM;
    qxM = qx(nodeInM);  qxP = 2.*hP.*phys.uin - qxM;
    qyM = qy(nodeInM);  qyP = 2.*qyM - qyM;
    [Fh, ~, ~] = SWE_Mex_BC2d...
        (hmin, gra, hM, hP, qxM, qxP, qyM, qyP, nx, ny);
    % assignment to numerical flux
    Fhs(mesh.mapI)  = Fh;
    Fqxs(mesh.mapW) = Fqx;
    Fqys(mesh.mapW) = Fqy;
end% if

if (strncmp(phys.casename, 'TsuamiRunup', 11))
    if (time<phys.inWave(1,end))
        nx  = mesh.nx(mesh.mapI);
        ny  = mesh.ny(mesh.mapI);
        % input water height
        hin = interp1(phys.inWave(1,:),phys.inWave(2,:),time,'linear');
        hM  = h(nodeInM);   hP  = hin - bot(nodeInM);
        qxM = qx(nodeInM);  qxP = qxM;
        qyM = qy(nodeInM);  qyP = 0;

        [Fh, ~, ~] = SWE_Mex_BC2d...
            (hmin, gra, hM, hP, qxM, qxP, qyM, qyP, nx, ny);
        % assignment to numerical flux
        Fhs(mesh.mapI)  = Fh;
    end% if
    return;
end% if

end% func