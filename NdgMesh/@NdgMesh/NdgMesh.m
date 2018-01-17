%> @brief Unstructed mesh object
%
%> NdgMesh class stores the information of the unstructed mesh,
%> including the elements and interpolation nodes connection.
%> Each NdgMesh objects only contain one kind of std cell. For the hybrid
%> mesh, the user have to setup multiple NdgMesh objects and connect them
%> with the NdgEdge object by
%> @code Matlab
%>   mesh1 = NdgMesh2d(tri, Nvt, vxt, vyt, vzt, Kt, EToVt, EToRt, BCToVt); 
%>   mesh2 = NdgMesh2d(quad, Nvq, vxq, vyq, vzq, Kq, EToVq, EToRq, BCToVq); 
%>   % first to connect the mesh object 
%>   mesh1.assembleMeshConnection(mesh2, 1, 2);
%>   mesh2.assembleMeshConnection(mesh1, 2, 1);
%>   % create the edge connection in the mesh object 
%>   mesh1.assembleNdgEdgeConnection(mesh2, 1, 2);
%>   mesh2.assembleNdgEdgeConnection(mesh1, 2, 1);
%> @endcode
%> 
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgMesh < handle
    
    properties
        %> std cell object
        cell
        %> number of cell
        K
        %> number of vertices
        Nv
        %> vertex index in each cell (column)
        EToV
        %> region types for each cell
        EToR int8
        %> coordinate of vertex
        vx
        %> coordinate of vertex
        vy
        %> coordinate of vertex
        vz
        %> edge objects
        edgeUnion
    end
    
    % elemental volume infomation
    properties( SetAccess=protected )
        %> boundary types for each cell (column)
        EToB
        %> mesh id of adjacent cell
        EToM
        %> adjacent cell index for each cell
        EToE
        %> adjacent face index for each cell
        EToF
        %> coordinate of interpolation points
        x
        %> coordinate of interpolation points
        y
        %> coordinate of interpolation points
        z
        %> determination of Jacobian matrix at each interpolation points
        J
        %> 
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        %> length/area/volume of each cell
        LAV
        %> character length of each cell
        charLength
    end
    
    properties%( SetAccess = protected )
        %> the local face point index
        eidM
        %> the adjacent point index
        eidP
        %> edge type of each facial points
        eidtype int8
    end
    
    properties( SetAccess = protected )
        %> central coordinate
        xc
        %> central coordinate
        yc
        %> central coordinate
        zc
        %> normal vector of each facial point
        nx
        %> normal vector of each facial point
        ny
        %> normal vector of each facial point
        nz
        %> determination of facial integral at each face points
        Js
    end
    
    properties( Hidden = true, SetAccess = protected )
        %> figure handle
        figureHandle  
    end
    
    methods( Abstract, Hidden, Access = protected )
        %> Get volume infomation of each element
%         [ J ] = assembleJacobiFactor( obj )
        [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobiFactor( obj )
        %> Get outward normal vector of each elemental edges
        [ nx, ny, nz ] = assembleFacialJaobiFactor( obj )
        [ faceId ] = assembleGlobalFaceIndex( obj )
        [ edge ] = makeConnectNdgEdge( obj, mesh1, mid0, mid1 )
    end
    
    methods(Hidden, Access=protected)
        [ EToE, EToF, EToM ] = assembleCellConnect( obj )
        [ eidM, eidP, eidtype ] = assembleEdgeNode( obj )
        [ EToB ] = assembleCellBoundary( obj, BCToV )
        
        function [LAV, charLength] = assembleCellScale( obj, J )
            
            % Jacobian determination on each quadrature points
            one = ones( obj.cell.Np, obj.K );
            LAV = mxGetMeshIntegralValue( one, obj.cell.wq, J, obj.cell.Vq );
            %LAV = mxGetIntegralValue(, obj.cell.wq, obj.Jq );
            switch obj.cell.type
                case NdgCellType.Line
                    charLength = LAV;
                case NdgCellType.Tri
                    charLength = sqrt( 2*LAV );
                case NdgCellType.Quad
                    charLength = sqrt( LAV );
            end
        end% func
        
        function [x, y, z] = assembleNodeCoor( obj, vx, vy, vz )
            x = obj.proj_vert2node(vx);
            y = obj.proj_vert2node(vy);
            z = obj.proj_vert2node(vz);
        end% func
    end
    
    methods( Abstract )
        %> determine the cell index that the gauge points in.
        [ cellId ] = accessGaugePointLocation( obj, xg, yg, zg );
    end
    
    methods
        
        function obj = NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV)
            [ obj.cell, obj.Nv, obj.vx, obj.vy, obj.vz, obj.K, obj.EToV, obj.EToR ] ...
                = checkInput(cell, Nv, vx, vy, vz, K, EToV, EToR);
            
            [ obj.EToE, obj.EToF, obj.EToM ]  = assembleCellConnect( obj );
            [ obj.x, obj.y, obj.z ] = assembleNodeCoor( obj, vx, vy, vz );
            [ obj.EToB ] = assembleCellBoundary(obj, BCToV);
            
            [ obj.rx, obj.ry, obj.rz, ...
                obj.sx, obj.sy, obj.sz, ...
                obj.tx, obj.ty, obj.tz, obj.J ] = assembleJacobiFactor( obj );
%             [ obj.J ] = assembleJacobiFactor( obj );
            [ obj.LAV, obj.charLength ] = assembleCellScale( obj, obj.J );
            [ obj.nx, obj.ny, obj.nz, obj.Js ] = assembleFacialJaobiFactor( obj );
            [ obj.eidM, obj.eidP, obj.eidtype ] = assembleEdgeNode( obj );
            
            [ obj.xc ] =  obj.GetMeshAverageValue( obj.x );
            [ obj.yc ] =  obj.GetMeshAverageValue( obj.y );
            [ obj.zc ] =  obj.GetMeshAverageValue( obj.z );
        end% func
        
        function assembleNdgEdgeConnection( obj, mesh1, mid0, mid1 )
            edge = obj.makeConnectNdgEdge( mesh1, mid0, mid1 );
            if edge.M > 0
                obj.edgeUnion = [ obj.edgeUnion, edge ];
            end
        end
        
        function assembleMeshConnection( obj, mesh1, mid0, mid1 )
            obj.EToM(:) = mid0;
            Nface = obj.cell.Nface;
            for n = 1:( obj.K * Nface )
                [f, k] = ind2sub( [Nface, obj.K], n );
                if ( obj.EToE(f, k) ~= k ) continue; end
                vert = sort( obj.EToV( obj.cell.FToV(:, f), k ) );
                
                Nface1 = mesh1.cell.Nface;
                for n1 = 1:( mesh1.K * Nface1 )
                    [f1, k1] = ind2sub( [Nface1, mesh1.K], n1 );
                    vert1 = sort( mesh1.EToV( mesh1.cell.FToV(:, f1), k1 ) );
                    if vert == vert1
                        obj.EToM(f, k) = mid1;
                        obj.EToE(f, k) = k1;
                        obj.EToF(f, k) = f1;
                        obj.setEToB(k, f, NdgEdgeType.GaussEdge)
                        break;
                    end
                end
            end
        end% func
        
        function nodeQ = proj_vert2node(obj, vertQ)
            % project scalars from mesh verts to nodes
            ele_vQ = vertQ(obj.EToV);
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
        
        function integralValue = GetMeshIntegralValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue...
                (nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
        end
        
        function avergeValue = GetMeshAverageValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue...
                (nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
            avergeValue = integralValue./obj.LAV;
        end
        
        function setEToB(obj, k, f, edgeType)
            obj.EToB(f, k) = NdgEdgeType( edgeType );
            endNfp = sum( obj.cell.Nfp(1:f) );
            facePointId = endNfp:-1:(endNfp - obj.cell.Nfp(f) + 1);
            obj.eidtype(facePointId, k) = int8( edgeType );
        end% func
        
    end% methods
    
end

function [cell, Nv, vx, vy, vz, K, EToV, EToR] = checkInput(cell, Nv, vx, vy, vz, K, EToV, EToR)

if( numel(vx) ~= Nv ) || (numel(vy) ~= Nv ) || (numel(vz) ~= Nv )
    msgID = [mfilename, ':InputVertexError'];
    msgtext = 'The number of input vertex is not equal to Nv.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if isrow(vx)
    vx = vx';
end

if isrow(vy)
    vy = vy';
end

if isrow(vz)
    vz = vz';
end

if( size(EToV, 1) ~= cell.Nv )
    msgID = [mfilename, ':InputCellError'];
    msgtext = 'The rows of EToV is not equal to Nv (cell).';
    ME = MException(msgID, msgtext);
    throw(ME);
end

if( size(EToV, 2) ~= K ) || (numel(EToR) ~= K)
    msgID = [mfilename, ':InputCellError'];
    msgtext = 'The number of cells in EToV or EToR is not equal to K.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

EToR = int8( NdgRegionType( EToR ) );
end
