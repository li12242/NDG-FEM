function [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, mesh1, mesh2, mid1, mid2 )

Nedge = 0;
ind = zeros(2, mesh1.K*mesh1.cell.Nface);

for n1 = 1:(mesh1.K*mesh1.cell.Nface)
    [f1, k1] = ind2sub([mesh1.cell.Nface, mesh1.K], n1);
    
    mshId = mesh1.EToM( n1 );
    if ( mshId ~= mid2 ) continue; end
    vert1 = sort( getFaceVertex(mesh1, k1, f1) );
    
    for n2 = 1:(mesh2.K*mesh2.cell.Nface)
        [f2, k2] = ind2sub([mesh2.cell.Nface, mesh2.K], n2);
        mshId = mesh2.EToM( n2 );
        if ( mshId ~= mid1 ) continue; end
        vert2 = sort( getFaceVertex(mesh2, k2, f2) );
        if ( vert1 == vert2 )
            Nedge = Nedge + 1;
            ind(1, Nedge) = n1;
            ind(2, Nedge) = n2;
        end
    end
end

FToE = zeros(2, Nedge);
FToF = zeros(2, Nedge);
FToV = zeros( max( mesh1.cell.Nfv ), Nedge );
% edgeType = NdgEdgeType.Inner * ones( Nedge, 1, 'int8' );
for n = 1:Nedge
    n1 = ind(1, n);
    n2 = ind(2, n);
    
    [f1, k1] = ind2sub([mesh1.cell.Nface, mesh1.K], n1);
    vert1 = sort( getFaceVertex(mesh1, k1, f1) );
    [f2, k2] = ind2sub([mesh2.cell.Nface, mesh2.K], n2);
    FToE(:, n) = [k1, k2]';
    FToF(:, n) = [f1, f2]';
    FToV(:, n) = vert1;
end
end% func

function vert = getFaceVertex( mesh, k, f )
locVertId = mesh.cell.FToV(:, f);
vert = mesh.EToV(locVertId, k);
end