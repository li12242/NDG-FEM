function ulimit = BJ2D(mesh, u, eleIndex)
% Barth and Jepson slope limiter
% Input:
%   mesh     - mesh object
%   u        - original variables
%   eleIndex - element index
% Output:
%   ulimit   - limited result
% 

V = mesh.Shape.VandMatrix;
Np = mesh.Shape.nNode;

%% Compute cell averages
% The cell average is obtained by dividing the Vandermonde matrix
uh         = V\u;  
uh(2:Np,:) = 0; 
uavg       = V*uh;  
v          = uavg(1,:);

%% Compute limited result

% Initialize the limited result
ulimit = u;

% find max and min cell averages
maxv = max(v(mesh.EToE(eleIndex, :)'));
minv = min(v(mesh.EToE(eleIndex, :)'));

maxv = max([v(eleIndex); maxv]);
minv = min([v(eleIndex); minv]);

% Apply reconstruction to elements
maxu  = ones(Np, 1)*maxv; 
minu  = ones(Np, 1)*minv;
meanu = ones(Np, 1)*v(:, eleIndex);
  
% correction factor
a = limitCoeff(maxu, minu, meanu, u(:, eleIndex));

% apply slope limiter to selected elements
ulimit(:,eleIndex) = meanu + a.*(ulimit(:, eleIndex) - meanu);
end% func


function a = limitCoeff(maxu, minu, meanu, u)
Np = size(u, 1);

a = min(1, ( maxu - meanu )./( u - meanu ) );

ind = u < meanu;
a(ind) = min(1, ( minu(ind) - meanu(ind) )./( u(ind) - meanu(ind) ) );

amin = min(a);
a = ones(Np, 1)*amin;
end