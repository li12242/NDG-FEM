%% PositivePreserving
% The operator to preserve the positive of wet elements.
% Refer to Xing (2010) and Zhang (2010) for more details.
function [h, qx, qy] = SWE_PositivePreserving2d(phys, mesh, h, qx, qy)

% Parameters
minDepth  = phys.minDepth;
% ksi       = 0.0;
% Np        = mesh.Shape.nNode;
% 
% % Compute mean water depth
% hmean     = CellMean(mesh, h);
% 
% % Correct mean water less than hDelta
% dis       = (hmean <= ksi);
% h(:, dis) = h(:, dis) + ones(Np, 1)*(ksi - hmean(dis));
% qx(:, dis)= 0;
% qy(:, dis)= 0;
% 
% hmean     = CellMean(mesh, h);
% qxmean    = CellMean(mesh, qx);
% qymean    = CellMean(mesh, qy);
% 
% % Obrain positive operator
% hmin      = min(h);
% theta     = min( (hmean - ksi)./(hmean - hmin), 1);
% theta(theta<0) = 0; % in case for hmean = hmin
% 
% % Reconstruction
% h  = (ones(Np, 1)*theta).*(h - ones(Np, 1)*hmean) + ones(Np, 1)*hmean;
% qx = (ones(Np, 1)*theta).*(qx - ones(Np, 1)*qxmean) + ones(Np, 1)*qxmean;
% qy = (ones(Np, 1)*theta).*(qy - ones(Np, 1)*qymean) + ones(Np, 1)*qymean;
% 
% flag = h<=minDepth;
% qx(flag) = 0;
% qy(flag) = 0;
% 
% % eliminate dry flux
% hmean     = CellMean(mesh, h);
% dis       = hmean<minDepth;
% qx(:,dis) = 0;
% qy(:,dis) = 0;

[h, qx, qy] = SWE_Mex_PositivePreserving2d(h, qx, qy, mesh.Shape.M, ...
    mesh.J, minDepth);
end% func

%% Compute cell averages
% The cell average is obtained by dividing the Vandermonde matrix.
function v = CellMean(mesh, u)
V  = mesh.Shape.V;
Np = mesh.Shape.nNode;

uh         = V\u;  
uh(2:Np,:) = 0; 
uavg       = V*uh;  
v          = uavg(1,:);
end