function disFlag = DisDetector(mesh, var, u)
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

disFlag = zeros(size(mesh.nx));

% get inflow boundary
inflowFlag = (mesh.nx.*uM < 0 );

disFlag(inflowFlag) = abs(vM(inflowFlag) - vP(inflowFlag));

% mean value
varMean = mesh.Shape.VandMatrix\var;
varMean = [varMean(1,:); varMean(1,:)];

% mean value
varFlag = (varMean <= eps);

h = mesh.x(end, :) - mesh.x(1, :);
h = [h; h]; % element length
% h = ones(size(mesh.nx));

disFlag = disFlag./(h.^((mesh.Shape.nOrder+1)/2).*varMean );
disFlag(varFlag) = 0;
disFlag = sum(disFlag, 1);
end% func