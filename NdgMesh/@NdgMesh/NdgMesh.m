%> @brief Unstructed mesh object
% ======================================================================
%> NdgMesh class stores the information of unstructed mesh,
%> including the elements and interpolation nodes connection.
%> Each NdgMesh objects only contain one kind of ref cell. For the hybrid
%> mesh, the user have to setup multiple NdgMesh objects and connect them
%> with other objects by
%> @code Matlab
%>   mesh1 = NdgMesh2d(tri, Nvt, vxt, vyt, vzt, Kt, EToVt, EToRt, BCToVt);
%>   mesh2 = NdgMesh2d(quad, Nvq, vxq, vyq, vzq, Kq, EToVq, EToRq, BCToVq);
%>   % first to connect the mesh object
%>   mesh1.assembleMeshConnection(mesh2);
%>   mesh2.assembleMeshConnection(mesh1);
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
        %> vertex index in each cell
        EToV
        %> region id for each cell
        EToR
        %> coordinate of vertex
        vx, vy, vz
        %> edge objects
        InnerEdge % inner edge
        BoundaryEdge % halo edge
        %> mesh index
        ind
        %> element status
        status      int8
    end
    
    % elemental volume infomation
    properties ( SetAccess=protected )
        %> mesh id of adjacent cell
        EToM
        %> adjacent cell index for each cell
        EToE
        %> adjacent face index for each cell
        EToF
        %> coordinate of interpolation points
        x, y, z
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
    
    properties ( SetAccess = protected )
        %> central coordinate
        xc, yc, zc
        %> normal vector of each facial point
        %nx, ny, nz
        %> determination of facial integral at each face points
        %Js
    end
    
    methods ( Abstract, Hidden, Access = protected )
        %> Get volume infomation of each element
        obj = assembleJacobiFactor( obj )
        %> Get outward normal vector of each elemental edges
        % [ nx, ny, nz ] = assembleFacialJaobiFactor( obj )
        % [ faceId ] = assembleGlobalFaceIndex( obj )
    end
    
    methods ( Hidden, Access=protected )
        obj = ConnectMesh( obj )
        % [ eidM, eidP, eidtype ] = assembleEdgeNode( obj )
        % [ EToB ] = assembleCellBoundary( obj, BCToV )
        obj = GetNodeCoor( obj )
        obj = GetCellSize( obj )
    end
    
    methods
        
        function obj = NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR)
            obj = checkInput( obj, cell, Nv, vx, vy, vz, K, EToV, EToR );
            obj = ConnectMesh( obj );
            obj = GetNodeCoor( obj );
            
            obj = assembleJacobiFactor( obj );
            obj = GetCellSize( obj );
            
            [ obj.xc ] =  obj.GetMeshAverageValue( obj.x );
            [ obj.yc ] =  obj.GetMeshAverageValue( obj.y );
            [ obj.zc ] =  obj.GetMeshAverageValue( obj.z );
        end% func
        
        % assemble mesh connection
        ConnectMeshUnion( obj, meshId, meshUnion )

        % project scalars from mesh verts to nodes
        function nodeQ = proj_vert2node(obj, vertQ)
            ele_vQ = vertQ(obj.EToV);
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
        
        function integralValue = GetMeshIntegralValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue(nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
        end
        
        function avergeValue = GetMeshAverageValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue(nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
            avergeValue = integralValue ./ obj.LAV;
        end
        
        [ err ] = evaluateNormErr1( obj, fphys, fext );
        [ err ] = evaluateNormErr2( obj, fphys, fext );
        [ err ] = evaluateNormErrInf( obj, fphys, fext );
    end% methods
end

function obj = checkInput(obj, cell, Nv, vx, vy, vz, K, EToV, EToR)

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

% assignment
obj.cell = cell;
obj.Nv = Nv;
obj.vx = vx;
obj.vy = vy;
obj.vz = vz;
obj.K = K;
obj.EToV = EToV;
obj.EToR = EToR;

end
