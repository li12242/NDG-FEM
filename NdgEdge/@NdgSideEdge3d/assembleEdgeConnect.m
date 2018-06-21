function obj = assembleEdgeConnect( obj, mesh )

edge2d = mesh.mesh2d.InnerEdge;
Nz = mesh.Nz;
Ne = edge2d.Ne * Nz;

obj.Ne = Ne;
% connect element
obj.FToE = [];
for n = 1 : Nz
    ind1 = ( edge2d.FToE(1,:) - 1) * Nz + n;
    ind2 = ( edge2d.FToE(2,:) - 1) * Nz + n;
    obj.FToE = [ obj.FToE; ind1 ];
    obj.FToE = [ obj.FToE; ind2 ];
end
obj.FToE = reshape( obj.FToE, 2, Ne );

obj.FToF = repmat( edge2d.FToF, Nz, 1 );
obj.FToF = reshape( obj.FToF, 2, Ne );

obj.FToM = mesh.ind;
obj.FToV = zeros( 4, obj.Ne );
for n = 1 : Ne
    e1 = obj.FToE( 1, n );
    f1 = obj.FToF( 1, n );
    i1 = mesh.cell.FToV(:, f1);
    obj.FToV( :, n ) = mesh.EToV( i1, e1 );
end
end% func