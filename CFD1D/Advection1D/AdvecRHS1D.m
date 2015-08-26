function [rhsu] = AdvecRHS1D(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

% Globals1D;

% form field differences at faces
alpha=0;
Mesh = u.mesh; Element = Mesh.Element;

% du = zeros(Mesh.Element.Nfp*Mesh.Element.Nfaces, Mesh.K);
du = u.evaluateJump.*(a*Mesh.nx(:)-(1-alpha)*abs(a*Mesh.nx(:)))/2;

% impose boundary condition at x=0
uin = -sin(a*time);
du (Mesh.mapI) = (u.val(Mesh.vmapI)- uin).*(a*Mesh.nx(Mesh.mapI)-(1-alpha)*...
    abs(a*Mesh.nx(Mesh.mapI)))/2;
du (Mesh.mapO) = 0;
du = reshape(du, Mesh.Element.Nfp*Mesh.Element.Nfaces, Mesh.K);

% compute right hand sides of the semi-discrete PDE
rhsu = -a*Mesh.rx.*(Element.Dr*u.val) + Element.LIFT*(Mesh.Fscale.*(du));
return
