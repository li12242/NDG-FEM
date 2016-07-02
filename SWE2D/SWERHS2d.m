function [rhsH, rhsQx, rhsQy] = SWERHS2d(phys, mesh, h, qx, qy)
% Calculation the RHS of SWE

% parameters
minDepth = phys.minDepth;
shape    = mesh.Shape;
bot      = phys.bot;
% find the dry nodes
dryEleFlag  = IsDry(mesh, h, minDepth);
dryNodeFlag = repmat(dryEleFlag, mesh.Shape.nNode, 1);

% flux
[Fh, Fqx, Fqy, Gh, Gqx, Gqy] = SWEFlux2d(phys, h, qx, qy, dryNodeFlag);
% source terms
[Sh, Sqx, Sqy] = SWESource2d(phys, mesh, h, bot, dryEleFlag);

% numerical flux
[Fhs, Fqxs, Fqys] = SWENumFlux2d(phys, mesh, h, qx, qy, dryEleFlag);

% the flux difference
dFh   = Fh(mesh.vmapM).*mesh.nx + Gh(mesh.vmapM).*mesh.ny   - Fhs;
dFqx  = Fqx(mesh.vmapM).*mesh.nx + Gqx(mesh.vmapM).*mesh.ny - Fqxs;
dFqy  = Fqy(mesh.vmapM).*mesh.nx + Gqy(mesh.vmapM).*mesh.ny - Fqys;

% dFqx(dFqx<1e-10) = 0;
% dFqy(dFqx<1e-10) = 0;

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