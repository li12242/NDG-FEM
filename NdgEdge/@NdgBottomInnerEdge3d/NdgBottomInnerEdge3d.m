classdef NdgBottomInnerEdge3d < handle
    %NDGHORIZONEDGE3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( SetAccess = protected )
        %> mesh obj
        mesh
        %> num of face nodes
        Nfp
        %> num of edges
        Ne
        %> mass matrix of edge
        M
        %> vertex index on each edge
        FToV
        %> local and adjacent cell index
        FToE
        %> local face index of local and adjacent cell
        FToF
        %> face to mesh index
        FToM
        %> interp node index of 1st ele on each edge
        FToN1
        %> interp node index of 2nd ele on each edge
        FToN2
        %> outward normal vector
        nx, ny, nz
        %> determination of edge Jacabian
        Js
%         %> edge type
%         ftypes
    end
    
    methods (Access = public)
        function obj = NdgBottomInnerEdge3d( meshUnion3d, meshId )
            mesh = meshUnion3d( meshId );

            obj.mesh = mesh;
            obj = assembleMassMatrix( obj, mesh.cell.N, mesh.cell.Nz );
            obj = assembleEdgeConnect( obj, mesh, mesh.mesh2d );
            obj = assembleNodeProject( obj, mesh );
        end

        %> access boundary values at edges
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
        %> evaluate strong-form surface term rhs
        [ frhs ] = matEvaluateStrongFormEdgeRHS( obj, fluxM, fluxP, fluxS )
        [ frhs ] = matEvaluateStrongFormEdgeAlterRHS( obj, fluxM, fluxP )
        [ frhs ] = matEvaluateStrongFormEdgeCentralRHS( obj, fluxM, fluxP )
    end
    
    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh, mesh2d );
        obj = assembleNodeProject( obj, mesh );
        obj = assembleMassMatrix( obj, N, Nz );
    end
end

