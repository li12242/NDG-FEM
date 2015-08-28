function [rhsH, rhsQ] = SWERHS(mesh, h, q, bedElva)
% 1D shallow water equation 
% Righ Hand Side

line = mesh.Shape;

% numel flux
[Fhs, Fqs] = SWEHLL(mesh, h, q);

% boundary condition
% [Fhs, Fqs] = SWEBC(mesh, Fhs, Fqs, h, q);

% source term
[Sh, Sq] = SWESource(mesh, bedElva, h);

% flux term
[Fh, Fq] = SWEFlux(h, q);

%% strong form
dFh = mesh.nx.*Fh(mesh.vmapM) - Fhs;
dFq = mesh.nx.*Fq(mesh.vmapM) - Fqs;
% 
% % rhs
rhsH = -mesh.rx.*(line.Dr*Fh) + line.invM*line.Mes*(dFh.*mesh.fScale) + Sh;
rhsQ = -mesh.rx.*(line.Dr*Fq) + line.invM*line.Mes*(dFq.*mesh.fScale) + Sq;

%% weak form
% rhs
% rhsH = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fh))) ...
%     - line.invM*line.Mes*(Fhs.*mesh.fScale) + Sh;
% rhsQ = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fq))) ...
%     - line.invM*line.Mes*(Fqs.*mesh.fScale) + Sq;

% eliminate the dry area
dryPointFlag = (h <= 10^-6);
eleFlag = (dryPointFlag(mesh.vmapM) > 0) & (dryPointFlag(mesh.vmapP) >0);
eledist = find(all(eleFlag));

rhsH(:,eledist) = zeros(mesh.Shape.nNode, numel(eledist));
rhsQ(:,eledist) = zeros(mesh.Shape.nNode, numel(eledist));
end% func


function [Fh, Fq] = SWEFlux(h, q)
g = 9.8; hDelta = 10^-6;
u = zeros(size(h)); flag = h>hDelta;
u(flag) = q(flag)./h(flag);
Fh = q;
Fq = g*h.^2./2 + u.^2.*h;
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
% % Inflow
% q0 = 0.18;
% Fhs(mesh.mapI) = q0.*mesh.nx(mesh.mapI);
% % Outflow
% h0 = 0.33; u0 = q(vmapO)./h0;
% Fqs(mesh.mapO) = (g*h0.^2./2 + u0.^2.*h0).*mesh.nx(mesh.mapO);

end% func

function [Sh, Sq] = SWESource(mesh, bedElva, h)
g = 9.8; line = mesh.Shape;
Sh = zeros(size(h)); %Sq = zeros(size(h));
Sq = -g.*h.*mesh.rx.*(line.Dr*bedElva);
end% func

function [Fhs, Fqs] = SWEHLL(mesh, h, q)
% numerical flux: HLL
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 60-60. (4.46)
% [2] [Song_2010] An unstructured finite volume model for dam-break 
%     floods with wet/dry fronts over complex topography. (17-18)
g = 9.8; TOL=10e-5; hDelta = 10^-6;

% wave speed
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = q(mesh.vmapM); qP = q(mesh.vmapP);
uM = zeros(size(hM)); uP = zeros(size(hP));
flagM = hM>hDelta; flagP = hP>hDelta;
uM(flagM) = qM(flagM)./hM(flagM); 
uP(flagP) = qP(flagP)./hP(flagP);

cs = 0.5*( sqrt(g*hM) + sqrt(g*hP) ) + 0.25*(uM - uP).*mesh.nx;
us = ( 0.5*(uM + uP) ).*mesh.nx + (sqrt(g*hM) - sqrt(g*hP));

us(cs<0) = 0; % middle is dry 
cs(cs<0) = 0; % middle is dry

SL = mesh.nx.*(uM) - sqrt(g.*hM);
SR = mesh.nx.*(uP) + sqrt(g.*hP);

SL = min(SL, us - cs); SR = max(SR, us + cs);

% left is dry
isDryLeft = (hM<TOL);
TL = mesh.nx.*(uP) - 2*sqrt(g*hP);
SL(isDryLeft) = TL(isDryLeft);

% right is dry
isDryRight = (hP<TOL);
TR = mesh.nx.*(uM) + 2*sqrt(g*hM);
SR(isDryRight) = TR(isDryRight);

% 
[FhM, FqM] = SWEFlux(hM, qM);
[FhP, FqP] = SWEFlux(hP, qP);

% Compute HLL flux
SDelta = SR-SL;
SDelta(abs(SDelta)<10^-16) = 10^-16;

t1 = (min(SR,0)-min(0,SL))./SDelta;
t2 = 1-t1;
t3 = (SR.*abs(SL)-SL.*abs(SR))./(2*SDelta);

Fhs = t1.*(FhP.*mesh.nx) + t2.*(FhM.*mesh.nx) - t3.*(hP-hM);
Fqs = t1.*(FqP.*mesh.nx) + t2.*(FqM.*mesh.nx) - t3.*(qP-qM);
end% func