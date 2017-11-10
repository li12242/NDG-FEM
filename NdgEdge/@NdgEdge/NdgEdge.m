%> @brief Edge information on unstructured meshes.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgEdge < handle
    
    properties(SetAccess=protected)
        %> edge std cell for mesh1
        bcell1
        %> edge std cell for mesh2
        bcell2
        %> number of edges
        M
        %> vertex index on each edge
        FToV
        %> local and adjacent cell index
        FToE
        %> local face index of local and adjacent cell
        FToF
        %> adjacent mesh index
        FToM
        
        %> element node index of the 1st mesh on each edge
        FToN1
        %> element node index of the 2nd mesh on each edge
        FToN2
        
        
        %> project matrix from the 2nd mesh edge nodes to the 1st mesh edge nodes
        IntM
%         %> project the boundary nodal values to the quadrature points for adjacent cell
%         NToQ2
%         
%         %> determination of Jacobian matrixs
%         Jq
%         %> x component of the outward normal vector
%         nxq
%         %> y component of the outward normal vector
%         nyq
%         %> z component of the outward normal vector
%         nzq
    end
    
    methods(Abstract, Hidden=true, Access=protected)
        [ bcell1, bcell2 ] = setStdEdgeCell( obj, mesh1, mesh2, meshId1, meshId2 );
%         [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale( obj, meshArray );
    end
    
    methods
        function obj = NdgEdge( mesh1, mesh2, mid1, mid2 )
            [ obj.FToM ] = mid2;
            [ obj.bcell1, obj.bcell2 ] = setStdEdgeCell( obj, mesh1, mesh2 );
            [ obj.M, obj.FToE, obj.FToF, obj.FToV ] = assembleEdgeConnect( obj, mesh1, mesh2, mid1, mid2 );
            [ obj.IntM, obj.FToN1, obj.FToN2 ] = assembleNodeProject( obj, mesh1, mesh2 );
%             [ obj.nxq, obj.nyq, obj.nzq, obj.Jq ] = assembleQuadPointScale( obj, mesh1 );
        end
    end
    
    methods(Hidden, Access=protected)
        [ Nedge, FToE, FToF, FToV ]  = assembleEdgeConnect( obj, mesh1, mesh2, meshId1, meshId2 );
        [ intM, FToN1, FToN2 ] = assembleNodeProject( obj, mesh1, mesh2 )
    end
end

