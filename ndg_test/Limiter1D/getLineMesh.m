function mesh = getLineMesh(x1, x2, eleNum, nOrder)
% Input:
%   x1 - start point
%   x2 - end point
%   eleNum - element number
%   nOrder - order of accuracy
% Output:
%   mesh - mesh object


% max order of polymomials
N = nOrder; nElement = eleNum;
[Nv, VX, K, EToV] = Utilities.Mesh.MeshGen1D(x1, x2, nElement);
BC = [2,1; 3,Nv];

% Initialize solver and construct grid and metric
line = StdRegions.Line(N);
mesh = MultiRegions.RegionLineBC(line, EToV, VX, BC);

end% func