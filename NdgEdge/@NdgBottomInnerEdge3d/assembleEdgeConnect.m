function obj = assembleEdgeConnect( obj, mesh, mesh2d )

Nz = mesh.Nz;
Kh = mesh2d.K;
% from top to bottom layer
Ne = Kh * (Nz - 1);

obj.Ne = Ne;
% connect element
obj.FToE = [];
for n = 1 : (Nz - 1)
    % index of up layer element
    ind1 = ( (1:Kh) - 1) * Nz + n;
    % index of bottom layer element
    ind2 = ( (1:Kh) - 1) * Nz + n + 1;
    temp = [ ind1; ind2 ];
    obj.FToE = [ obj.FToE, temp ];
end

% last layer
% ind1 = ( (1:Kh) - 1) * Nz + Nz;
% temp = [ ind1; ind1 ];
% obj.FToE = [ obj.FToE, temp ];

% % all edge are default bottom edge
obj.FToF = ones( 2, Ne );
% find middle edge
% ind = ( obj.FToE(1, :) == obj.FToE(2, :) );
% obj.FToF(1, ~ind) = 4;
% obj.FToF(2, ~ind) = 5;
if mesh.cell.type == enumStdCell.PrismTri
    obj.FToF(1, :) = 4;
    obj.FToF(2, :) = 5;
elseif mesh.cell.type == enumStdCell.PrismQuad
    obj.FToF(1, :) = 5;
    obj.FToF(2, :) = 6;
end

obj.FToM = mesh2d.ind;
obj.FToV = [];
for n = 1 : Nz
    v = mesh2d.EToV + mesh2d.Nv * n;
    obj.FToV = [v, obj.FToV];
end

% ftype = zeros(obj.Ne, 1);
% bottomLayerId = Kh * (Nz - 1) + ( 1 : Kh );
% ftype( bottomLayerId ) = enumBoundaryCondition.BottomBoundary;
% obj.ftype = enumBoundaryCondition( ftype );
end% func