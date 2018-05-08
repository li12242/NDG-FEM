adv = ConstAdvUniformMesh2d( 1, 5, NdgCellType.Quad );
mesh = adv.meshUnion;
edge = NdgInnerEdge2d( mesh, 1 );

xfM = zeros( edge.Nfp, edge.Ne );
yfM = zeros( edge.Nfp, edge.Ne );
xfP = zeros( edge.Nfp, edge.Ne );
yfP = zeros( edge.Nfp, edge.Ne );

for i = 1:edge.Ne
    k1 = edge.FToE(1, i);
    k2 = edge.FToE(2, i);
    
    fp1 = edge.FToFN1(:, i);
    fp2 = edge.FToFN2(:, i);
    xfM(:, i) = mesh.x( mesh.cell.Fmask( fp1 ), k1 );
    xfP(:, i) = mesh.x( mesh.cell.Fmask( fp2 ), k2 );
    yfM(:, i) = mesh.y( mesh.cell.Fmask( fp1 ), k1 );
    yfP(:, i) = mesh.y( mesh.cell.Fmask( fp2 ), k2 );
end

