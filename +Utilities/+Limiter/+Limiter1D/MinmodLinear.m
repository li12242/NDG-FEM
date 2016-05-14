<<<<<<< HEAD
=======
<<<<<<< HEAD:+Utilities/+Limiter/+Limiter1D/MinmodLinear.m
function ulimit = MinmodLinear(mesh, u)
% function ulimit = SlopeLimit1(u);
% Purpose: Apply minmod limiter to u, get linear result

V = @(x)mesh.Shape.VandMatrix;
Np = @(x)mesh.Shape.nNode;

% Compute modal coefficients
uh = V()\u;

% Extract linear polynomial, ul
ul = uh; ul(3:Np(),:) = 0; ul = V()*ul;

% Extract cell averages
uh(2:Np(),:)=0; uavg = V()*uh; v = uavg(1,:);

% Find cell averages in neighborhood of each element
vk = v;
vkm1 = [v(1),v(1:mesh.nElement-1)]; 
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)]; 

% Limit function in all cells
ulimit = Utilities.Limiter.Limiter1D.LimitLinear(ul,mesh.x,vkm1,vk,vkp1, mesh);
return
=======
>>>>>>> master
function ulimit = MinmodLinear(mesh, u)
% function ulimit = SlopeLimit1(u);
% Purpose: Apply minmod limiter to u, get linear result

V = @(x)mesh.Shape.VandMatrix;
Np = @(x)mesh.Shape.nNode;

% Compute modal coefficients
uh = V()\u;

% Extract linear polynomial, ul
ul = uh; ul(3:Np(),:) = 0; ul = V()*ul;

% Extract cell averages
uh(2:Np(),:)=0; uavg = V()*uh; v = uavg(1,:);

% Find cell averages in neighborhood of each element
vk = v;
vkm1 = [v(1),v(1:mesh.nElement-1)]; 
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)]; 

% Limit function in all cells
ulimit = Utilities.Limiter.Limiter1D.LimitLinear(ul,mesh.x,vkm1,vk,vkp1, mesh);
return
<<<<<<< HEAD
=======
>>>>>>> 4c67ef677518787526fc4d2de5e78498d098910c:+Utilities/+Limiter/+Limiter1D/MinmodLinear.m
>>>>>>> master
