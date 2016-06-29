function [rhsH, rhsQx, rhsQy] = SWERHS2d(phys, mesh, h, qx, qy, dryEleFlag)
% Calculation the RHS of SWE

% parameters
minDepth = phys.minDepth;
shape    = mesh.Shape;
bot      = phys.bot;
% find the dry nodes

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

% weak form
divFh  = Div2D(mesh, Fh, Gh);
% divFh(:, dryEleFlag) = 0.0;  % eliminate the dry element flux terms
rhsH   = -divFh + Sh + shape.LIFT*(dFh.*mesh.fScale);

divFh  = Div2D(mesh, Fqx, Gqx);
rhsQx  = -divFh + Sqx + shape.LIFT*(dFqx.*mesh.fScale);

divFh  = Div2D(mesh, Fqy, Gqy);
rhsQy  = -divFh + Sqy + shape.LIFT*(dFqy.*mesh.fScale);
end