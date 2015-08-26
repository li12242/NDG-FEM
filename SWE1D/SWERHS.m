function [rhsH, rhsQ] = SWERHS(mesh, h, q, bedElva)
% Idealized dam break problem of 1D shallow water equation 
% Righ Hand Side

line = mesh.Shape;
% numel flux
[Fhs, Fqs] = SWEHLL(mesh, h, q);

% boundary condition
[Fhs, Fqs] = SWEBC(mesh, Fhs, Fqs, h, q);

% flux term
[Fh, Fq] = SWEFlux(h, q);

dFh = mesh.nx.*Fh(mesh.vmapM) - Fhs;
dFq = mesh.nx.*Fq(mesh.vmapM) - Fqs;

% source term
[Sh, Sq] = SWESource(mesh, bedElva, h);

% rhs
rhsH = -mesh.rx.*(line.Dr*Fh) + line.invM*line.Mes*(dFh.*mesh.fScale) + Sh;
rhsQ = -mesh.rx.*(line.Dr*Fq) + line.invM*line.Mes*(dFq.*mesh.fScale) + Sq;
end% func

function [Fh, Fq] = SWEFlux(h, q)
g = 9.8;
Fh = q;
Fq = g*h.^2./2 + q.^2./h;
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
g = 9.8; TOL=10e-5;

% wave speed
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = q(mesh.vmapM); qP = q(mesh.vmapP);

cs = 0.5*( sqrt(g*hM) + sqrt(g*hP) ) + 0.25*(qM./hM - qP./hP).*mesh.nx;
us = ( 0.5*(qM./hM + qP./hP) ).*mesh.nx + (sqrt(g*hM) - sqrt(g*hP));

us(cs<0) = 0; % middle is dry 
cs(cs<0) = 0; % middle is dry

SL = mesh.nx.*(qM./hM) - sqrt(g.*hM);
SR = mesh.nx.*(qP./hP) + sqrt(g.*hP);

SL = min(SL, us - cs); SR = max(SR, us + cs);

% left is dry
isDryLeft = (hM<TOL);
TL = mesh.nx.*(qP./hP) - 2*sqrt(g*hP);
SL(isDryLeft) = TL(isDryLeft);

% right is dry
isDryRight = (hP<TOL);
TR = mesh.nx.*(qM./hM) + 2*sqrt(g*hM);
SR(isDryRight) = TR(isDryRight);

% 
[FhM, FqM] = SWEFlux(hM, qM);
[FhP, FqP] = SWEFlux(hP, qP);

% Compute HLL flux
t1 = (min(SR,0)-min(0,SL))./(SR-SL);
t2 = 1-t1;
t3 = (SR.*abs(SL)-SL.*abs(SR))./(2*(SR-SL));

Fhs = t1.*(FhP.*mesh.nx) + t2.*(FhM.*mesh.nx) - t3.*(hP-hM);
Fqs = t1.*(FqP.*mesh.nx) + t2.*(FqM.*mesh.nx) - t3.*(qP-qM);
end% func