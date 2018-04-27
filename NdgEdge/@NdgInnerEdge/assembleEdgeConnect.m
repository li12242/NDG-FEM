
function [ Nedge, FToE, FToF, FToV, ftype ] = assembleEdgeConnect( obj, mesh1 )

ind = zeros(mesh1.K*mesh1.cell.Nface, 1);
for k = 1:mesh1.K
    for f = 1:mesh1.cell.Nface
        ind( f+(k-1)*mesh1.cell.Nface ) = getUniqueFaceIndex( mesh1, k, f );
    end
end

uniqueInd = unique( ind );
Nedge = numel( uniqueInd );

FToE = zeros(2, Nedge);
FToF = zeros(2, Nedge);
FToV = zeros( max( mesh1.cell.Nfv ), Nedge );
ftype = NdgEdgeType.Inner * ones( Nedge, 1, 'int8' );
for n = 1:Nedge
    id = uniqueInd( n );
    tid = find( ind == id );
    
    if ( numel(tid) == 2 )
        [f1, k1] = ind2sub([mesh1.cell.Nface, mesh1.K], tid(1));
        vert1 = sort( getFaceVertex(mesh1, k1, f1) );
        [f2, k2] = ind2sub([mesh1.cell.Nface, mesh1.K], tid(2));
        FToE(:, n) = [k1, k2]';
        FToF(:, n) = [f1, f2]';
        FToV(:, n) = vert1;
        ftype( n ) = mesh1.EToB(f1, k1);
    else % Boundary edge or Halo edge
        [f1, k1] = ind2sub([mesh1.cell.Nface, mesh1.K], tid );
        vert1 = sort( getFaceVertex(mesh1, k1, f1) );
        FToE(:, n) = [k1, k1]';
        FToF(:, n) = [f1, f1]';
        FToV(:, n) = vert1;
        type = mesh1.EToB(f1, k1);
        if ( type == NdgEdgeType.Inner ) % connect to other mesh elements
            ftype( n ) = NdgEdgeType.GaussEdge;
        else
            ftype( n ) = type;
        end
    end
end
end% func

function vert = getFaceVertex( mesh, k, f )
locVertId = mesh.cell.FToV(:, f);
vert = mesh.EToV(locVertId, k);
end

function ind = getUniqueFaceIndex( mesh, k, f )
locVertId = mesh.cell.FToV(:, f);
vert = sort( mesh.EToV(locVertId, k) );
ind = sum( vert .* (mesh.Nv .^ ((numel(vert)-1):-1:0)') );
end
