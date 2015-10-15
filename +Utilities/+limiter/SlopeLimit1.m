function ulimit = SlopeLimit1(mesh, u)
% function ulimit = SlopeLimit1(u);
% Purpose: Apply slopelimiter (Pi^1) to u
V = @(x)mesh.Shape.VandMatrix;
Np = @(x)mesh.Shape.nNode;

% Compute modal coefficients
uh = V()\u;

% Extract linear polynomial
ul = uh; ul(3:Np(),:) = 0;ul = V()*ul;

% Extract cell averages
uh(2:Np(),:)=0; uavg = V()*uh; v = uavg(1,:);

% Find cell averages in neighborhood of each element
vk = v; vkm1 = [v(1),v(1:mesh.nElement-1)]; 
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)]; 

% Limit function in all cells
ulimit = Utilities.Limiter.SlopeLimitLin(ul,mesh.x,vkm1,vk,vkp1, mesh);
return
