function [rhsu] = PoissonCstabRHS1D(u, Element, Mesh)

% function [rhsu] = PoissonCstabRHS1D(u,time)
% Purpose  : Evaluate RHS in 1D Poisson equation on
%            symmetric form using stabilized central flux

% Globals1D;

% Define field differences at faces
Nfp = Element.Nfp; Nfaces = Element.Nfaces;
K = Mesh.K;

du = zeros(Nfp*Nfaces,K); du(:) = (u(Mesh.vmapM)-u(Mesh.vmapP))/2.0;

% impose boundary condition -- Dirichlet BC's
uin  = -u(Mesh.vmapI); du (Mesh.mapI) = (u(Mesh.vmapI) -  uin )/2.0;
uout = -u(Mesh.vmapO); du (Mesh.mapO) = (u(Mesh.vmapO) - uout)/2.0;

% Compute q
q = Mesh.rx.*(Element.Dr*u) - Element.LIFT*(Mesh.Fscale.*(Mesh.nx.*du));
dq = zeros(Nfp*Nfaces,K); dq(:) = q(Mesh.vmapM)-q(Mesh.vmapP);

% impose boundary condition -- Neumann BC's
qin  = q(Mesh.vmapI); dq (Mesh.mapI) = q(Mesh.vmapI)- qin; 
qout = q(Mesh.vmapO); dq (Mesh.mapO) = q(Mesh.vmapO)-qout;

% evaluate fluxes
tau = 1.0;
fluxq = Mesh.nx.*(dq/2.0+tau*Mesh.nx.*du);

% compute right hand sides of the semi-discrete PDE
rhsu = Mesh.J.*((Element.invV'*Element.invV)*(Mesh.rx.*(Element.Dr*q)...
    - Element.LIFT*(Mesh.Fscale.*fluxq)));
return
