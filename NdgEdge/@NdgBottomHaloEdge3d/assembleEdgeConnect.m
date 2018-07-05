function obj = assembleEdgeConnect( obj, mesh, mesh2d )
%ASSEMBLEEDGECONNECT Summary of this function goes here
%   Detailed explanation goes here

Kh = mesh2d.K;
Nz = mesh.Nz;
Ne = Kh;

obj.Ne = Ne;
% connect element

% last layer
ind1 = ( (1:Kh) - 1) * Nz + Nz;
obj.FToE  = [ ind1; ind1 ];

obj.FToF = ones( 2, Ne );
if mesh.cell.type == enumStdCell.PrismTri
    obj.FToF(1:2, :) = 4;
elseif mesh.cell.type == enumStdCell.PrismQuad
    obj.FToF(1:2, :) = 5;
end

obj.FToM = mesh2d.ind;
obj.FToV = mesh2d.EToV + mesh2d.Nv * Nz;

ftype = zeros(obj.Ne, 1);
ftype(:) = enumBoundaryCondition.BottomBoundary;
obj.ftype = enumBoundaryCondition( ftype );
end
