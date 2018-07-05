function obj = assembleEdgeConnect( obj, mesh )
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description
    
edge2d = mesh.mesh2d.BoundaryEdge;
Nz = mesh.Nz; % num of vertical layers
Ne = edge2d.Ne * Nz; % num of edge

obj.Ne = Ne;

% connect element
obj.FToE = [];
for n = 1 : Nz
    ind1 = ( edge2d.FToE(1,:) - 1) * Nz + n; % local cell index
    ind2 = ( edge2d.FToE(2,:) - 1) * Nz + n; % adjacent cell index
    obj.FToE = [ obj.FToE; ind1 ];
    obj.FToE = [ obj.FToE; ind2 ];
end
obj.FToE = reshape( obj.FToE, 2, Ne );

obj.FToF = reshape( repmat( edge2d.FToF, Nz, 1 ), 2, Ne);
obj.FToM = reshape( repmat( edge2d.FToM, Nz, 1 ), 2, Ne);

% edge is quadrilateral with 4 vertices
obj.FToV = zeros( 4, obj.Ne );
for n = 1 : Ne
    e1 = obj.FToE( 1, n );
    f1 = obj.FToF( 1, n );
    i1 = mesh.cell.FToV(:, f1);
    obj.FToV( :, n ) = mesh.EToV( i1, e1 );
end

obj.ftype = repmat( edge2d.ftype', Nz, 1 );
obj.ftype = obj.ftype(:);
end