function [Fhs, Fqs] = SWELF(mesh, h, q)
% Lax-Friedrichs numerical flux

hDelta = 10^-3; g = 9.8;

% rotate variable
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = mesh.nx.* q(mesh.vmapM); qP = mesh.nx.* q(mesh.vmapP);
uM = zeros(size(hM)); uP = zeros(size(hP));
isWetM = hM > hDelta; isWetP = hP > hDelta;

uM(isWetM) = qM(isWetM)./hM(isWetM);
uP(isWetP) = qP(isWetP)./hP(isWetP);

% estimate wave speed
SM = abs(uM) + sqrt(g*hM); SP = abs(uP) + sqrt(g*hP);
SM(~isWetM) = 0; SM(~isWetP) = 0;
S = max(SM, SP);

[FhM, FqM] = SWEFlux(hM, qM, hDelta);
[FhP, FqP] = SWEFlux(hP, qP, hDelta);

% Compute HLL flux
% Fhs = zeros(size(FhM)); Fqs = zeros(size(FqM));
Fhs = 0.5*(FhM + FhP + S.*(hM - hP));
Fqs = 0.5*(FqM + FqP + S.*(qM - qP));

% both dry
sflag = ~isWetM & ~isWetP;
Fhs(sflag) = 0; Fqs(sflag) = 0;

% rotate invariance
Fqs = Fqs.*mesh.nx;
end% func