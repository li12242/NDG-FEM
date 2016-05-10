function ulimit = MinmodModified(mesh, u)
% function ulimit = MinmodModified(u);
% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th order polynomial


V = @(x)mesh.Shape.VandMatrix;
Np = @(x)mesh.Shape.nNode;

% Compute cell averages, v
uh = V()\u; uh(2:Np(),:)=0; uavg = V()*uh; v = uavg(1,:);

% Apply slope limiter as needed
ulimit = u; eps0=1.0e-8;

% find end values of each element
ue1 = u(1,:); ue2 = u(end,:);

% find cell averages
vk = v;
vkm1 = [v(1),v(1:mesh.nElement-1)];
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)];

% Apply reconstruction to find elements in need of limiting
ve1 = vk - Utilities.Limiter.Limiter1D.MinmodFun([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + Utilities.Limiter.Limiter1D.MinmodFun([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

% Check to see if any elements require limiting
if(~isempty(ids))
  % create piecewise linear solution for limiting on specified elements
  uhl = V()\u(:,ids); uhl(3:Np(),:)=0; ul = V()*uhl;
  
  % apply slope limiter to selected elements
  ulimit(:,ids) = Utilities.Limiter.Limiter1D.LimitLinear(ul, ...
      mesh.x(:,ids),vkm1(ids),vk(ids),vkp1(ids),mesh);
end% if
end% func
