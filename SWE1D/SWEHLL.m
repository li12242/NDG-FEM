function [Fhs, Fqs, status] = SWEHLL(mesh, h, q)
% numerical flux: HLL
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 60-60. (4.46)
% [2] [Song_2010] An unstructured finite volume model for dam-break 
%     floods with wet/dry fronts over complex topography. (17-18)
% [3] 
% 
hDelta = 10^-3;

% rotate variable
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = mesh.nx.* q(mesh.vmapM); qP = mesh.nx.* q(mesh.vmapP);
uM = zeros(size(hM)); uP = zeros(size(hP));
isWetM = hM > hDelta; isWetP = hP > hDelta;

uM(isWetM) = qM(isWetM)./hM(isWetM);
uP(isWetP) = qP(isWetP)./hP(isWetP);

[SM, SP, status] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);

[FhM, FqM] = SWEFlux(hM, qM, hDelta);
[FhP, FqP] = SWEFlux(hP, qP, hDelta);

% Compute HLL flux
Fhs = zeros(size(FhM)); Fqs = zeros(size(FqM));

sflag = (SM >= 0);
Fhs(sflag) = FhM(sflag);
Fqs(sflag) = FqM(sflag);

sflag = SP <= 0;
Fhs(sflag) = FhP(sflag);
Fqs(sflag) = FqP(sflag);

sflag = SM < 0 & SP > 0;
Fhs(sflag) = (- SM(sflag).*FhP(sflag) ...
    + SP(sflag).*FhM(sflag)...
    + SP(sflag).*SM(sflag).*( hP(sflag) - hM(sflag) ) )...
    ./(SP(sflag) - SM(sflag));
Fqs(sflag) = (- SM(sflag).*FqP(sflag) ...
    + SP(sflag).*FqM(sflag)...
    + SP(sflag).*SM(sflag).*( qP(sflag) - qM(sflag) ) )...
    ./(SP(sflag) - SM(sflag));

% both dry
sflag = ~isWetM & ~isWetP;
Fhs(sflag) = 0; Fqs(sflag) = 0;

% rotate invariance
Fqs = Fqs.*mesh.nx;

end% func