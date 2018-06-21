classdef NdgSideEdge3d < handle
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
    end

    methods ( Access = public )
        function obj = NdgSideEdge3d( meshUnion3d, meshId )
            mesh = meshUnion3d( meshId );
            
            obj.mesh = mesh;
            obj = assembleMassMatrix( obj, mesh.cell.N, mesh.cell.Nz );
            obj = assembleEdgeConnect( obj, mesh );
            obj = assembleNodeProject( obj, mesh );
        end

        %> access boundary values at edges
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
        %> evaluate strong-form surface term rhs
        [ frhs ] = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS )
    end

    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
        [ Nfp, M ] = assembleMassMatrix( obj, N, Nz );

        [ nx, ny, nz, Js ] = PrismQuadJacobian3d( obj, mesh, f1, e1, fid );
        [ nx, ny, nz, Js ] = PrismTriJacobian3d( obj, mesh, f1, e1, fid );
    end
end