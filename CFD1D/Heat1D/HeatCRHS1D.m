function [rhsu] = HeatCRHS1D(u,time,Element,Mesh)

% function [rhsu] = HeatCRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D heat equation using central flux

% Globals1D;

% Define field differences at faces
du = zeros(Element.Nfp*Element.Nfaces,Mesh.K); 
du(:)  = (u(Mesh.vmapM)-u(Mesh.vmapP))/2.0;

% impose boundary condition -- Dirichlet BC's
uin  = -u(Mesh.vmapI); du(Mesh.mapI) = (u(Mesh.vmapI)-uin)/2.0; 
uout = -u(Mesh.vmapO); du(Mesh.mapO)=(u(Mesh.vmapO) - uout)/2.0;

% Compute q and form differences at faces
q = Mesh.rx.*(Element.Dr*u) - Element.LIFT*(Mesh.Fscale.*(Mesh.nx.*du));
dq = zeros(Element.Nfp*Element.Nfaces,Mesh.K); 
dq(:)  = (q(Mesh.vmapM)-q(Mesh.vmapP))/2.0;

% impose boundary condition -- Neumann BC's
qin  = q(Mesh.vmapI); dq(Mesh.mapI) = (q(Mesh.vmapI)- qin )/2.0; 
qout = q(Mesh.vmapO); dq(Mesh.mapO) = (q(Mesh.vmapO)-qout)/2.0;

% compute right hand sides of the semi-discrete PDE
rhsu = Mesh.rx.*(Element.Dr*q) - Element.LIFT*(Mesh.Fscale.*(Mesh.nx.*dq));
return
