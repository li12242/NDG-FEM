%======================================================================
%> @brief assemble mesh connection
%> @detals
%>
%> @param mesh1 other mesh
%>
%======================================================================
%> This function is part of the NDGOM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
%======================================================================
function assembleMeshConnection( obj, mesh1 )

% obj.EToM(:) = obj.ind; % local mesh index
% Nface = obj.cell.Nface;
% for n = 1:( obj.K * Nface )
%     [f, k] = ind2sub( [Nface, obj.K], n );
%     if ( obj.EToE(f, k) ~= k )
%         continue;
%     end
%     vert = sort( obj.EToV( obj.cell.FToV(:, f), k ) );
%     
%     Nface1 = mesh1.cell.Nface;
%     for n1 = 1:( mesh1.K * Nface1 )
%         [f1, k1] = ind2sub( [Nface1, mesh1.K], n1 );
%         vert1 = sort( mesh1.EToV( mesh1.cell.FToV(:, f1), k1 ) );
%         if vert == vert1
%             obj.EToM(f, k) = mesh1.ind;
%             obj.EToE(f, k) = k1;
%             obj.EToF(f, k) = f1;
%             obj.setEToB(k, f, NdgEdgeType.GaussEdge);
%             break;
%         end
%     end
% end

[ obj.EToM, obj.EToE, obj.EToF, obj.EToB ] ...
    = mxAssembleMeshConnection( obj.ind, mesh1.ind, ...
    obj.K, mesh1.K, obj.cell.Nface, mesh1.cell.Nface, ...
    obj.cell.FToV, mesh1.cell.FToV, ...
    obj.EToV, mesh1.EToV, obj.EToM, obj.EToE, obj.EToF, obj.EToB);

end% func