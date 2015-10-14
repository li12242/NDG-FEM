function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1,mesh)

% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right

% Globals1D;

Np = mesh.Shape.nNode;
% Compute various geometric measures
ulimit = ul; h = xl(Np,:)-xl(1,:); 
x0 = ones(Np,1)*(xl(1,:) + h/2);

hN = ones(Np,1)*h;

% Limit function
ux = (2./hN).*(mesh.Shape.Dr*ul);

ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*Utilities.minmod([ux(1,:);(vp1-v0)./h;(v0-vm1)./h]));
return
