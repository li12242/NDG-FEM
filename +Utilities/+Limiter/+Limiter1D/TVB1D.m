function ulimit = TVB1D(mesh, u, M)
% TVB slope limiter from Cockburn & Shu (1989)
% Input:
%   M - parameters

V = mesh.Shape.VandMatrix;
Np = mesh.Shape.nNode;

% Compute cell averages, v
uh = V\u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

% Apply slope limiter as needed
ulimit = u;

% find end values of each element
ue1 = u(1,:); ue2 = u(end,:);

% find cell averages
vk = v;
vkm1 = [v(1),v(1:mesh.nElement-1)];
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)];

% the threshold
l = (mesh.J(1, :)*2).^2*M; % Mh2

% Apply reconstruction to find elements in need of limiting
ids = find( abs(vk-ue1) > l | abs(ue2-vk) > l );

% Check to see if any elements require limiting
if(~isempty(ids))
  % create piecewise linear solution for limiting on specified elements
  uhl = V\u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
  
  % apply slope limiter to selected elements
  ulimit(:,ids) = Utilities.Limiter.Limiter1D.LimitLinear(ul, ...
      mesh.x(:,ids),vkm1(ids),vk(ids),vkp1(ids),mesh);
end% if

end% func