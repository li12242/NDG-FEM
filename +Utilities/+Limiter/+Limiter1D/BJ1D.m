function ulimit = BJ1D(mesh, u)

V = mesh.Shape.VandMatrix;
Np = mesh.Shape.nNode;

% Compute cell averages, v
uh = V\u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

% Apply slope limiter as needed
ulimit = u;

% find max and min cell averages
maxv = max(v(mesh.EToE')); 
minv = min(v(mesh.EToE')); 

maxv = max([v; maxv]);
minv = min([v; minv]);

% Apply reconstruction to find elements in need of limiting
minu = min(u);
maxu = max(u);

dx = mesh.J(1, :)*2; K = 100;

ids = find( (maxu > maxv | minu < minv) & (maxu - minu) > K.*dx);

% Check to see if any elements require limiting
if(~isempty(ids))
  
  % 
  maxu = ones(Np, 1)*maxv; minu = ones(Np, 1)*minv;
  meanu = ones(Np, 1)*v;
  
  % correction factor
  a = limitCoeff(maxu, minu, meanu, u);
  
  % apply slope limiter to selected elements
  ulimit(:,ids) = meanu(:, ids) + ...
      a(:, ids).*(ulimit(:, ids) - meanu(:, ids));
end% if

end% func


function a = limitCoeff(maxu, minu, meanu, u)
a = ones(size(u)); Np = size(u, 1);

ind = u > meanu;
a(ind) = min(1, ( maxu(ind) - meanu(ind) )./( u(ind) - meanu(ind) ) );

ind = u < meanu;
a(ind) = min(1, ( minu(ind) - meanu(ind) )./( u(ind) - meanu(ind) ) );

amin = min(a);
a = ones(Np, 1)*amin;
end