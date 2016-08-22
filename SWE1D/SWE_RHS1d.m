%% SWE_RHS1d
% Calculate the right hand side of one dimensional SWE.
function [rhsH, rhsQ] = SWE_RHS1d(phys, mesh, h, q, bot, isWet, time, dt)
% Righ Hand Side
line = mesh.Shape;
% numel flux
[Fhs, Fqs]  = SWE_HLL1d(phys, mesh, h, q, isWet);
% source term
[Sh, Sq]    = SWE_Source1d(phys, mesh, bot, h, isWet);
% flux term
[Fh, Fq]    = SWE_Flux1d(phys, h, q, isWet);

[Sph, Spq]  = SWE_SpongeLayerBC(phys, mesh, h, q, time, dt);

% rhs
rhsH = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fh))) ...
    - line.LIFT*(Fhs.*mesh.fScale) + Sh + Sph;
rhsQ = mesh.rx.*(line.invM*(line.Dr'*(line.M*Fq))) ...
    - line.LIFT*(Fqs.*mesh.fScale) + Sq + Spq;
end% func

%% SWE_SpongeLayerBC
% set external solution condition for sponge layer.
function [Sph, Spq]  = SWE_SpongeLayerBC(phys, mesh, h, q, time, dt) 
if isempty(phys.spl)
    Sph = zeros(size(mesh.x));
    Spq = zeros(size(mesh.x));
else
    for i = 1:numel(phys.spl)
        sigma = phys.spl(i).GetAbsorpCoeff(mesh, dt);
        temp  = phys.spl(i).GetBC(mesh, 'hb', time);
        Sph   = -sigma.*(h - temp);
        temp  = phys.spl(i).GetBC(mesh, 'qb', time);
        Spq   = -sigma.*(q - temp);
    end
end% if
end% func

%% SWE_LF1d
% Lax Lax-Friedrichs numerical flux
function [Fhs, Fqs] = SWE_LF1d(phys, mesh, h, q, isWet)

hDelta = phys.minDepth; 
g = 9.8;
% rotate variable
hM = h(mesh.vmapM); hP = h(mesh.vmapP);
qM = mesh.nx.* q(mesh.vmapM); qP = mesh.nx.* q(mesh.vmapP);
uM = zeros(size(hM)); uP = zeros(size(hP));
isWetM = hM > hDelta; isWetP = hP > hDelta;

uM(isWetM) = qM(isWetM)./hM(isWetM);
uP(isWetP) = qP(isWetP)./hP(isWetP);

% estimate wave speed
SM = abs(uM) + sqrt(g*hM); SP = abs(uP) + sqrt(g*hP);
SM(~isWetM) = 0; 
SM(~isWetP) = 0;
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

%% SWE_HLL1d
% Calculation of HLL function.
function [Fhs, Fqs] = SWE_HLL1d(phys, mesh, h, q, isWet)
% REFERENCE:
% [1] Khan A A, Lai W. Modeling Shallow Water Flows Using the 
%     Discontinuous Galerkin Method[M]. CRC Press, 2014. 60-60. (4.46)
% [2] [Song_2010] An unstructured finite volume model for dam-break 
%     floods with wet/dry fronts over complex topography. (17-18)
% 
% mex version
hmin = phys.minDepth;
gra  = phys.gra;
[Fhs, Fqs] = SWE_Mex_HLLFlux1d(hmin, gra, h, q, mesh.nx, mesh.vmapM, mesh.vmapP);


% minDepth = phys.minDepth;
% 
% % rotate variable
% hM      = h(mesh.vmapM); 
% hP      = h(mesh.vmapP);
% qM      = mesh.nx.* q(mesh.vmapM); 
% qP      = mesh.nx.* q(mesh.vmapP);
% uM      = zeros(size(hM)); 
% uP      = zeros(size(hP));
% dryM    = hM > minDepth;
% dryP    = hP > minDepth;
% 
% uM(dryM) = qM(dryM)./hM(dryM);
% uP(dryP) = qP(dryP)./hP(dryP);
% [SM, SP] = SWE_HLLWaveSpeed1d(phys, hM, hP, uM, uP);
% 
% % [FhM, FqM] = SWEFlux(hM, qM, isWet, physics);
% [FhM, FqM] = SWE_Flux1d(phys, hM, qM, isWet);
% % [FhP, FqP] = SWEFlux(hP, qP, isWet, physics);
% [FhP, FqP] = SWE_Flux1d(phys, hP, qP, isWet);
% 
% % Compute HLL flux
% Fhs = zeros(size(FhM)); 
% Fqs = zeros(size(FqM));
% 
% sflag = (SM >= 0); 
% Fhs(sflag) = FhM(sflag); Fqs(sflag) = FqM(sflag);
% 
% sflag = (SP <= 0);
% Fhs(sflag) = FhP(sflag); Fqs(sflag) = FqP(sflag);
% 
% sflag = ((SM < 0) & (SP > 0));
% Fhs(sflag) = (- SM(sflag).*FhP(sflag) + SP(sflag).*FhM(sflag)...
%     + SP(sflag).*SM(sflag).*(hP(sflag)-hM(sflag)))./(SP(sflag) - SM(sflag));
% Fqs(sflag) = (- SM(sflag).*FqP(sflag) + SP(sflag).*FqM(sflag)...
%     + SP(sflag).*SM(sflag).*( qP(sflag) - qM(sflag) ) )...
%     ./(SP(sflag) - SM(sflag));
% 
% % both dry
% sflag = ((~dryM) & (~dryP));
% Fhs(sflag) = 0; Fqs(sflag) = 0;
% 
% % rotate invariance
% Fqs = Fqs.*mesh.nx;
end% func

%% SWE_HLLWaveSpeed1d
% Wave speed of HLL funcion
function [SM, SP] = SWE_HLLWaveSpeed1d(phys, hM, hP, uM, uP)
% Estimate wave speed
gra      = phys.gra;
minDepth = phys.minDepth;
    
us       = 0.5*(uM + uP) + (sqrt(gra*hM) - sqrt(gra*hP));
cs       = 0.5*( sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(uM - uP);
% allocate mem
SM       = zeros(size(uM)); 
SP       = zeros(size(uP));
wetM     = (hM >= minDepth); 
wetP     = (hP >= minDepth);
% left is dry
flag     =  wetP & ~wetM; 
SM(flag) = uP(flag) - 2*sqrt(gra.*hP(flag) ); % hM = 0
SP(flag) = uP(flag) + sqrt(gra.*hP(flag) ); % hM = 0

% right is dry
flag     = ~wetP &  wetM; 
SP(flag) = uM(flag) + 2*sqrt(gra.*hM(flag) ); % hP = 0
SM(flag) = uM(flag) - sqrt(gra.*hM(flag) ); % hP = 0

% both wet
flag     =  wetP & wetM; 
SM(flag) = uM(flag) - sqrt(gra.*hM(flag));
SM(flag) = min(SM(flag), us(flag) - cs(flag) );
SP(flag) = uP(flag) + sqrt(gra.*hP(flag));
SP(flag) = max(SP(flag), us(flag) + cs(flag) );

end% func

%% SWE_Flux1d
% calculate flux terms of SWE
function [Fh, Fq] = SWE_Flux1d(phys, h, q, isWet)

hmin = phys.minDepth;
gra  = phys.gra;

[Fh, Fq] = SWE_Mex_Flux1d(hmin, gra, h, q);

% % Parameters
% minDepth = phys.minDepth;
% g        = phys.gra;
% % Wet/Dry status
% wetNode  = (h > minDepth);
% % flow rate
% u        = zeros(size(h)); 
% u(wetNode) = q(wetNode)./h(wetNode);
% 
% % Flux terms
% Fh = q;
% Fq = g*h.^2./2 + u.^2.*h;
% % for dry elements, no flow flux
% Fq(:, ~isWet) = 0.0;
end% func


%% SWE_Source1d
% Calculation of the source term.
function [Sh, Sq] = SWE_Source1d(phys, mesh, bot, h, isWet)

g    = phys.gra; 
line = mesh.Shape;
Sh   = zeros(size(h));
Sq   = -g.*h.*mesh.rx.*(line.Dr*bot);
Sq(:, ~isWet) = 0.0;
end% func



