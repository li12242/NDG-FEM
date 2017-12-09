classdef NdgGQMesh2d < NdgMesh2d
    
    methods(Hidden, Access=protected)
        function [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] ...
                = assembleJacobiFactor(obj)
            
            [rx, ry, rz, sx, sy, sz, tx, ty, tz, J] = ...
                assembleJacobiFactor@NdgMesh2d( obj );
            
            rx = obj.cell.project_node2quad( rx );
            ry = obj.cell.project_node2quad( ry );
            sx = obj.cell.project_node2quad( sx );
            sy = obj.cell.project_node2quad( sy );
            
            J = obj.cell.project_node2quad( J );
        end
        
        function [nx, ny, nz, Js] = assembleFacialJaobiFactor( obj )
            [nx, ny, nz, Js] = assembleFacialJaobiFactor@NdgMesh2d( obj );
            
            line = StdLine( obj.cell.N );
            for f = 1:obj.cell.Nface
                ind = sum( obj.cell.Nfp(1:f) );
                fpInd = (ind - obj.cell.Nfp(f) + 1):ind;
                nx(fpInd, :) = line.project_node2quad( nx(fpInd, :) );
                ny(fpInd, :) = line.project_node2quad( ny(fpInd, :) );
                nz(fpInd, :) = line.project_node2quad( nz(fpInd, :) );
                Js(fpInd, :) = line.project_node2quad( Js(fpInd, :) );
            end
        end
    end
    
    methods
        function obj = NdgGQMesh2d(cell, Nv, vx, vy, K, EToV, EToR, BCToV)
            obj = obj@NdgMesh2d(cell, Nv, vx, vy, K, EToV, EToR, BCToV);
        end% func
    end
    
end
