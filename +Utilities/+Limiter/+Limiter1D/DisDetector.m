function [disFlag, I] = DisDetector(mesh, var, u)
% Reference: 
%   [1]: Krivodonova 2004

% discontinuity detector
% Input:
%   mesh - mesh object
%   var - detecroted variable
%   u - flow speed 
% Output:
%   disFlag - bool flag for mesh with discontinuity

uM = u(mesh.vmapM); 
% get boundary value and neighbour's 
vM = var(mesh.vmapM); vP = var(mesh.vmapP);

I = zeros(size(mesh.nx));

% get inflow boundary
inflowFlag = (mesh.nx.*uM < 0 );

I(inflowFlag) = abs(vM(inflowFlag) - vP(inflowFlag));

% mean value
varMean = mesh.Shape.VandMatrix\var;
varMean = [varMean(1,:); varMean(1,:)];

% mean value
varFlag = ( abs(varMean) <= eps);

h = mesh.x(end, :) - mesh.x(1, :);
h = [h; h]; % element length
% h = ones(size(mesh.nx));

I = I./(h.^((mesh.Shape.nOrder+1)/2).*varMean );
I(varFlag) = 0;
I = sum(I, 1);

disFlag = I > 5;
end% func