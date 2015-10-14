function [Fhs, Fqs] = SWEHLL(mesh, h, q)
% numerical flux: HLL
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 60-60. (4.46)
% [2] [Song_2010] An unstructured finite volume model for dam-break 
%     floods with wet/dry fronts over complex topography. (17-18)
% [3] 
% 
hDelta = 10^-6;

% wave speed
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = q(mesh.vmapM); qP = q(mesh.vmapP);
uM = zeros(size(hM)); uP = zeros(size(hP));

isWetM = hM>hDelta; isWetP = hP>hDelta;

uM(isWetM) = qM(isWetM)./hM(isWetM);
uP(isWetP) = qP(isWetP)./hP(isWetP);

[SM, SP] = EstimateWaveSpeed(mesh, hM, hP, uM, uP);

[FhM, FqM] = SWEFlux(hM, qM);
[FhP, FqP] = SWEFlux(hP, qP);

% Compute HLL flux
Fhs = zeros(size(FhM)); Fqs = zeros(size(FqM));

sflag = SM >= 0;
Fhs(sflag) = FhM(sflag).*mesh.nx(sflag);
Fqs(sflag) = FqM(sflag).*mesh.nx(sflag);

sflag = SP <= 0;
Fhs(sflag) = FhP(sflag).*mesh.nx(sflag);
Fqs(sflag) = FqP(sflag).*mesh.nx(sflag);

sflag = SM < 0 & SP > 0;
Fhs(sflag) = (- SM(sflag).*(FhP(sflag).*mesh.nx(sflag)) ...
    + SP(sflag).*(FhM(sflag).*mesh.nx(sflag))...
    + SP(sflag).*SM(sflag).*( hP(sflag)-hM(sflag) ) )...
    ./(SP(sflag) - SM(sflag));
Fqs(sflag) = (- SM(sflag).*(FqP(sflag).*mesh.nx(sflag)) ...
    + SP(sflag).*(FqM(sflag).*mesh.nx(sflag))...
    + SP(sflag).*SM(sflag).*( qP(sflag) - qM(sflag) ) )...
    ./(SP(sflag) - SM(sflag));

sflag = ~isWetM & ~isWetP;
Fhs(sflag) = 0;
Fqs(sflag) = 0;
end% func