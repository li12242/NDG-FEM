function [rhsQ] = CNSRHS2D(Q, mu, time, SolutionBC, fluxtype)

% function [rhsQ] = CurvedCNSRHS2D(Q, mu, time, SolutionBC, fluxtype)
% Purpose: evaluate right hand side residual of the 
%          compressible Navier-Stokes equations
Globals2D;

% Gas constant
gamma = 1.4; 
vmapM = reshape(vmapM, Nfp*Nfaces, K); vmapP = reshape(vmapP, Nfp*Nfaces, K);

% ---------------------------------------------------------------------------

% Extract fields from three dimensional array
[~,~,rho,rhou,rhov,Ener] = CNSFlux2D(Q, gamma);

% ---------------------------------------------------------------------------

% Compute primitive fields at Gauss quadrature surface nodes
u = rhou./rho; v = rhov./rho; Pr = (gamma-1)*(Ener-0.5*rho.*(u.^2+v.^2));  

% ---------------------------------------------------------------------------

% Compute gradients of the conserved variables
[drhodx,   drhody] = Grad2D(rho);
[drhoudx, drhoudy] = Grad2D(rhou);
[drhovdx, drhovdy] = Grad2D(rhov);

% ---------------------------------------------------------------------------

% Use product-rule to evaluate gradients of velocity components at cubature nodes
dudx = (drhoudx  - drhodx.*u)./rho;
dudy = (drhoudy  - drhody.*u)./rho;
dvdx = (drhovdx  - drhodx.*v)./rho;
dvdy = (drhovdy  - drhody.*v)./rho;

% Compute viscous stress tensor at cubature nodes
ct11 = mu.*(2*dudx - (2/3)*(dudx + dvdy)); 
ct12 = mu.*(dudy + dvdx);
ct22 = mu.*(2*dvdy - (2/3)*(dudx + dvdy)); 
ct31 = u.*ct11 + v.*ct12; 
ct32 = u.*ct12 + v.*ct22; 

% ---------------------------------------------------------------------------

% Add mass conservation terms together and compute divergence
rhsQ(:,:,1) = - Div2D(rhou, rhov);

% Add x-momentum conservation terms together and compute divergence
rhsQ(:,:,2) = Div2D( - rhou.*u - Pr + ct11 , -rhou.*v + ct12);

% Add y-momentum conservation terms together and compute divergence
rhsQ(:,:,3) = Div2D( -rhou.*v + ct12, - rhov.*v - Pr + ct22);

% Add Energy conservation terms together and compute divergence
rhsQ(:,:,4) = Div2D(- u.*(Ener + Pr) + ct31, -v.*(Ener +Pr)+ct32);

% ---------------------------------------------------------------------------

% 2.1 evaluate '-' and '+' traces of conservative variables
for n=1:4
  Qn = Q(:,:,n);
  QM(:,:,n) = Qn(vmapM); QP(:,:,n) = Qn(vmapP);
end

% 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
[fM,gM,rhoM,rhouM,rhovM,EnerM] = CNSFlux2D(QM, gamma);
[fP,gP,rhoP,rhouP,rhovP,EnerP] = CNSFlux2D(QP, gamma);

% 2.2 set boundary conditions by modifying positive traces
if(~isempty(SolutionBC))
  [rhoT, rhouT, rhovT, EnerT] = feval(SolutionBC, rhoP, rhouP, rhovP, EnerP, time);
end

rhoP(mapB) = 2*rhoT(mapB) - rhoM(mapB); rhouP(mapB) = 2*rhouT(mapB) - rhouM(mapB); 
rhovP(mapB) = 2*rhovT(mapB) - rhovM(mapB); EnerP(mapB) = 2*EnerT(mapB) - EnerM(mapB); 

QP(:,:,1) = rhoP; QP(:,:,2) = rhouP; QP(:,:,3) = rhovP; QP(:,:,4) = EnerP;

% Add Lax-Friedrichs jump stabilization
lambda = sqrt(u.^2 + v.^2) + sqrt(abs(gamma*Pr./rho));
lambda = max(lambda(vmapM), lambda(vmapP));
lambda = reshape(lambda, Nfp, Nfaces*K);
lambda = ones(Nfp, 1)*max(lambda, [], 1);
lambda = reshape(lambda, Nfp*Nfaces, K);

for n=1:4
  nflux = (nx.*(fM(:,:,n) - fP(:,:,n)) + ny.*(gM(:,:,n) - gP(:,:,n)) - ...
      lambda.*(QM(:,:,n) - QP(:,:,n)))./2;
  rhsQ(:,:,n) = rhsQ(:,:,n) + LIFT*(Fscale.*nflux/2);
end

% ---------------------------------------------------------------------------
return;
