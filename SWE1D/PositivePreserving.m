function [h, q] = PositivePreserving(mesh, h, q, bedElva, isWet)
% Slope limiter and Positivity-preserving operator

%% slope limiter on water level and discharge

% the slope limiter act on the wet cells
q = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,q);

% eta = h + bedElva;
% eta = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,eta);

% temp = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,eta); 
% eta(:, wetIndex) = temp(:, wetIndex); % reconstruct dry element to linear

% h = eta - bedElva;

h = Utilities.Limiter.Limiter1D.MinmodLinear(mesh,h); 
% h(:, ~wetIndex) = temp(:, ~wetIndex);

%% positive preserving operator
[h, q] = PositiveOperator(mesh, h, q);

end% func

function [h, q] = PositiveOperator(mesh, h, q)
% positive operator
% reference from Xing (2010); Zhang (2010)

ksi = 0.0;

hmean = CellMean(mesh, h);
Np = mesh.Shape.nNode;
% correct mean water less than hDelta
dis = (hmean <= ksi);
h(:, dis) = h(:, dis) + ones(Np, 1)*(ksi - hmean(dis));
q(:, dis) = 0;

hmean = CellMean(mesh, h);
qmean = CellMean(mesh, q);
% positive operator
hmin = min(h);
theta = min( (hmean - ksi)./(hmean - hmin), 1);
h = (ones(Np, 1)*theta).*(h - ones(Np, 1)*hmean) + ones(Np, 1)*hmean;
q = (ones(Np, 1)*theta).*(q - ones(Np, 1)*qmean) + ones(Np, 1)*qmean;
end% func