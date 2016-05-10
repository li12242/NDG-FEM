function ulim = MMPlim(mesh, u)
% modification of the MP limiter

% Reference:
%   [1]: QIU J, SHU C W, 2005.

V = @(x)mesh.Shape.VandMatrix;
Np = @(x)mesh.Shape.nNode;
ulim = u;
eps0 = 1.0e-8;

% Compute modal coefficients
uh = V()\u;

% Extract cell averages
uh(2:Np(),:)=0; uavg = V()*uh; v = uavg(1,:);

vm1 = [v(1),v(1:mesh.nElement-1)]; 
vp1 = [v(2:mesh.nElement),v(mesh.nElement)];

dumin = v - min(min(vm1, v), vp1); 

% edge value
ue1 = u(1, :); ue2 = u(end, :);

dminu = v - min(ue1, ue2);


ids = find(abs(min(1, dumin./dminu) - 1) < eps0); % trouble cell flag

if(~isempty(ids))
  % create piecewise linear solution for limiting on specified elements
  uhl = V()\u(:,ids); uhl(3:Np(),:)=0; ul = V()*uhl;
  
  % apply slope limiter to selected elements
  ulim(:,ids) = Utilities.Limiter.Limiter1D.LimitLinear(ul, ...
      mesh.x(:,ids),vm1(ids),v(ids),vp1(ids),mesh);
end% if

end% function
