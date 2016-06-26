function ulimit = GKTVB1D(mesh, u)
% modified TVB limiter from Ghostine and Kesserwani (2009)

V = mesh.Shape.VandMatrix;
Np = mesh.Shape.nNode;

% Compute cell averages, v
uh = V\u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

% Apply slope limiter as needed
ulimit = u;

% find cell averages
vk = v;
vkm1 = [v(1),v(1:mesh.nElement-1)];
vkp1 = [v(2:mesh.nElement),v(mesh.nElement)];

% find max and min cell averages
maxv = max(v(mesh.EToE')); 
minv = min(v(mesh.EToE')); 

maxv = max([v; maxv]);
minv = min([v; minv]);

% Apply reconstruction to find elements in need of limiting
minu = min(u);
maxu = max(u);
ids = find( maxu > maxv | minu < minv );

% Check to see if any elements require limiting
if(~isempty(ids))
  % create piecewise linear solution for limiting on specified elements
  uhl = V\u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
  
  % apply slope limiter to selected elements
  ulimit(:,ids) = Utilities.Limiter.Limiter1D.LimitLinear(ul, ...
      mesh.x(:,ids),vkm1(ids),vk(ids),vkp1(ids),mesh);
end% if

end% func