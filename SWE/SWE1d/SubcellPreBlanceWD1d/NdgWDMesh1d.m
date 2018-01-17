classdef NdgWDMesh1d < NdgMesh1d
    
    properties( SetAccess = protected )
        Kmax
        coarseCellId
    end
    
    methods
        function obj = NdgWDMesh1d( mesh )
            cell = StdLine(1);
            obj = obj@NdgMesh1d( ...
                cell, mesh.Nv, mesh.vx, mesh.K, ...
                mesh.EToV, mesh.EToR, [] );
            obj.Kmax = mesh.K;
            obj.K = 0;
            obj.coarseCellId = zeros( mesh.K, 1 );
        end
        
        function reinitSubcellMesh( obj, meshId )
            obj.K = 0;
            obj.EToR(:) = int8( NdgRegionType.Dry );
            obj.EToM(:) = meshId;
        end
        
        function setSubcell( obj, vx, eidP, EToR, EToB, meshId, cellId )
            k = obj.K + 1;
            ele_vx = [vx(1), vx(3); vx(3), vx(2)];
            obj.x(:, [k, k+1]) = obj.cell.project_vert2node(ele_vx);
            obj.eidP(1, k  ) = eidP(1); % set adjacent node id
            obj.eidP(2, k+1) = eidP(2); 
            obj.eidtype(1, k  ) = EToB(1); % set element boundary node type
            obj.eidtype(2, k+1) = EToB(2); 
            obj.EToR( [k,k+1] ) = EToR(:);
            obj.EToM( [k,k+1] ) = meshId;
            obj.coarseCellId( [k, k+1] ) = cellId;
            
            xr = obj.cell.Dr*obj.x(:, [k, k+1]);
            obj.J(:, [k, k+1]) = xr;
            obj.Js(:, [k, k+1]) = 1;
            obj.rx(:, [k, k+1]) = 1./xr;
            [obj.LAV([k, k+1]), ~] = obj.assembleCellScale( obj.J(:, [k, k+1]) );
            obj.xc( [k, k+1] ) = obj.GetMeshIntegralValue( ...
                obj.J(:, [k, k+1]), ...
                obj.LAV([k, k+1]), ...
                obj.x(:, [k, k+1]) );
            obj.nx( :, [k, k+1] ) = sign( ...
                obj.x(obj.cell.Fmask, [k, k+1]) - obj.xc( [k, k+1] ) );
            % add element number
            obj.K = obj.K+2;
        end
        
        function [ fint ] = GetMeshIntegralValue(obj, J, LAV, node_val)
            fq = obj.cell.Vq * node_val;
            Jq = obj.cell.Vq * J;
            fint = (obj.cell.wq' * (Jq.*fq))./LAV;
        end
    end
    
    methods(Hidden, Access=protected)
        function [ LAV, charLength] = assembleCellScale( obj, J )
            Jq = obj.cell.Vq * J;
            LAV = obj.cell.wq' * Jq;
            charLength = LAV;
        end
    end
    
end
