function obj = assembleEdgeConnect( obj, mesh, mesh2d )

Nz = mesh.Nz;
Kh = mesh2d.K;
% from top to bottom layer
Ne = Kh * Nz;

obj.Ne = Ne;
% connect element
obj.FToE = [];
for n = 1 : (Nz - 1)
    % index of up layer element
    ind1 = ( (1:Kh) - 1) * Nz + n;
    % index of bottom layer element
    ind2 = ( (1:Kh) - 1) * Nz + n + 1;
    temp = [ind1; ind2];
    obj.FToE = [ obj.FToE, temp ];
end

% last layer
ind1 = ( (1:Kh) - 1) * Nz + Nz;
temp = [ind1; ind1];
obj.FToE = [ obj.FToE, temp ];

% all edge are default bottom edge
obj.FToF = ones( 2, Ne ) * 4;
% find middle edge
ind = ( obj.FToE(1, :) == obj.FToE(2, :) );
obj.FToF(1, ~ind) = 4;
obj.FToF(2, ~ind) = 5;

obj.FToM = mesh2d.ind;
obj.FToV = [];
for n = 1 : Nz
    v = mesh2d.EToV + mesh2d.Nv * n;
    obj.FToV = [v, obj.FToV];
end

end% func