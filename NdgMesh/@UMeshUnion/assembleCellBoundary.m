function [ EToB ] = assembleCellBoundary( obj, BCToV )
Nface = obj.cell.Nface;
Ke = obj.K;
Nbc = size(BCToV, 2);
EToB = NdgEdgeType.Inner * ones(obj.cell.Nface, obj.K, 'int8');
Nfp = size(BCToV, 1) - 1;

bcInd = false(Nbc, 1);
for k = 1:Ke
    for f = 1:Nface
        locVert = obj.cell.FToV(:, f);
        vert = obj.EToV(locVert, k);
        if ( obj.EToE(f, k) ~= k ) continue; end
        
        for b = 1:Nbc
            bvert = sort( BCToV(1:Nfp, b) );
        
            if ( sum( abs( bvert - sort(vert) ) ) < 1e-10 )
                EToB(f, k) = NdgEdgeType( BCToV(Nfp+1, b) );
                bcInd( b ) = true;
                break;
            end
        end
    end
end

if ( any( ~bcInd ) )
    bcInd = find( ~bcInd );
    msgID = 'UMeshUnion:assembleCellBoundary';
    msgtext = ['The edge of boundary condition #', num2str(bcInd'), ...
        ' is not found.'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func