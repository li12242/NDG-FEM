function [disFlag, I] = DisDetector(mesh, var, u, v)
% Reference: 
%   1. Krivodonova, Xin & Remacle et al. (2004)

% discontinuity detector
% Input:
%   mesh - mesh object
%   var - detecroted variable
%   u - flow speed 
% Output:
%   disFlag - bool flag for mesh with discontinuity

%% I. edge integral
shape = mesh.Shape;

% average matrix
AVE = sum(shape.Mes);

% outflow node index
ind = (mesh.nx.*u(mesh.vmapM) + mesh.ny.*v(mesh.vmapM) ) >= 0;

% eliminate outflow boundary value
temp = var(mesh.vmapM) - var(mesh.vmapP);
temp(ind) = 0;

% boundary integral of $\int_{\partial \Omega} l_{j} C_j$
I = AVE*(mesh.sJ.*abs(temp));

%% II. radius of the circumscribed circle
% Area of the standard element
area = sum(sum(shape.M));

r = sqrt(mesh.J(1,:)*area./pi);
I = I./r;%.^( (mesh.Shape.nOrder + 1)./2 );

%% III. edge length
temp = mesh.sJ; temp(ind) = 0;
l = AVE*temp;

% small = (abs(l) <= eps);
I = I./l;

%% IV. max norm
TOL = 1e-4;
% max norm
temp = max(abs(var));
ind = abs(temp) > TOL;

I(ind) = I(ind)./temp(ind);
I(~ind) = 0;

disFlag = I > 1;

end% func