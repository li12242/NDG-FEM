function [rhsH, rhsQ, status] = SWERHS(mesh, h, q, bedElva)
% 1D shallow water equation 
% Righ Hand Side
line = mesh.Shape;
hDelta = 1e-3;

% numel flux
% [Fhs, Fqs] = SWELF(mesh, h, q);
[Fhs, Fqs, status] = SWEHLL(mesh, h, q);


% boundary condition
% [Fhs, Fqs] = SWEBC(mesh, Fhs, Fqs, h, q);
% source term
[Sh, Sq] = SWESource(mesh, bedElva, h);

% flux term
[Fh, Fq] = SWEFlux(h, q, hDelta);

%% strong form
dFh = mesh.nx.*Fh(mesh.vmapM) - Fhs;
dFq = mesh.nx.*Fq(mesh.vmapM) - Fqs;

% eliminate dry boundary flux
isdry = (h <= hDelta);
dryEdge = (isdry(mesh.vmapM) & isdry(mesh.vmapP));
dFh(dryEdge) = 0; dFq(dryEdge) = 0;

% rhs
rhsH = -mesh.rx.*(line.Dr*Fh) + line.invM*line.Mes*(dFh.*mesh.fScale) + Sh;
rhsQ = -mesh.rx.*(line.Dr*Fq) + line.invM*line.Mes*(dFq.*mesh.fScale) + Sq;

%% weak form
% rhs
% rhsH = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fh))) ...
%     - line.invM*line.Mes*(Fhs.*mesh.fScale) + Sh;
% rhsQ = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fq))) ...
%     - line.invM*line.Mes*(Fqs.*mesh.fScale) + Sq;

%% eliminate the dry area gravitational flow
dryIndex = all(isdry);

rhsH(:,dryIndex) = line.invM*line.Mes*(dFh(:, dryIndex).*mesh.fScale(:, dryIndex));
rhsQ(:,dryIndex) = line.invM*line.Mes*(dFq(:, dryIndex).*mesh.fScale(:, dryIndex));

end% func

function [Sh, Sq] = SWESource(mesh, bedElva, h)
g = 9.8; line = mesh.Shape;
Sh = zeros(size(h)); %Sq = zeros(size(h));
Sq = -g.*h.*mesh.rx.*(line.Dr*bedElva);
end% func

function [Fhs, Fqs] = SWEBC(mesh, Fhs, Fqs, h, q)
g = 9.8; vmapI = mesh.vmapM(mesh.mapI); vmapO = mesh.vmapM(mesh.mapO);

%% subcritical flow
% % Inflow
% q0 = 0.18; 
% Fhs(mesh.mapI) = q0.*mesh.nx(mesh.mapI);
% 
% % Outflow
% h0 = 0.5; u0 = q(vmapO)./h0;
% Fqs(mesh.mapO) = (g*h0.^2./2 + u0.^2.*h0).*mesh.nx(mesh.mapO);

%% supercritical flow
% % Inflow
% q0 = 20.0567; h0 = 2.0;
% Fhs(mesh.mapI) = q0.*mesh.nx(mesh.mapI);
% Fqs(mesh.mapI) = (g*h0.^2./2 + q0.^2./h0).*mesh.nx(mesh.mapI);

%% transcritical flow
% Inflow
q0 = 0.18;
Fhs(mesh.mapI) = q0.*mesh.nx(mesh.mapI);
% Outflow
h0 = 0.33; u0 = q(vmapO)./h0;
Fqs(mesh.mapO) = (g*h0.^2./2 + u0.^2.*h0).*mesh.nx(mesh.mapO);

end% func



