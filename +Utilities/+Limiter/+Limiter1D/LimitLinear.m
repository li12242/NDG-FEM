function ulimit = LimitLinear(ul,xl,vm1,v0,vp1,mesh)
% function ulimit = LimitLinear(ul,xl,vm1,v0,vp1);
% Purpose: Apply minmod limited on linear function ul(Np,Ne) on x(Np,Ne)

% Input:
%   ul - linear polynomial
%   xl - node coordinate
%   vm1,v0,vp1 - cell averages left, center, and right
%   mesh - mesh object


Np = mesh.Shape.nNode; % number of nodes in an element

% Compute various geometric measures
h = xl(Np,:)-xl(1,:); 
x0 = ones(Np,1)*(xl(1,:) + h/2); % centre of element

hN = ones(Np,1)*h;

% gradient at each node, as ul is linear, each node has the same value
ux = (2./hN).*(mesh.Shape.Dr*ul); 

% Limit function
ulimit = ones(Np,1)*v0 + (xl-x0).*(ones(Np,1)...
    *Utilities.Limiter.Limiter1D.minmod([ux(1,:);(vp1-v0)./h;(v0-vm1)./h]));
return
