classdef NdgHaloEdge3d < NdgSideEdge3d
    properties ( SetAccess = protected )
        %> face type of each edge
        ftype
        %> boundary coordinate
        xb, yb, zb
    end

    methods (Access = public)
        function obj = NdgHaloEdge3d( meshUnion3d, meshId )
            obj = obj@NdgSideEdge3d(meshUnion3d, meshId);
        end

        %> evaluate right-hand-side for surface integral term
        frhs = matEvaluateStrongFormEdgeRHS( obj, fluxM, fluxP, fluxS );
        %> get surface values from physical field
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
    end
    
    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
    end

end