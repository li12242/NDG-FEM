classdef NdgBottomHaloEdge3d < NdgBottomInnerEdge3d
    %NDGBOTTOMHALOEDGE3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( SetAccess = protected )
        %> edge type
        ftype
    end
    
    methods
        function obj = NdgBottomHaloEdge3d( meshUnion3d, meshId )
            obj = obj@NdgBottomInnerEdge3d( meshUnion3d, meshId );
        end
        
        [ frhs ] = matEvaluateStrongFormEdgeRHS( obj, fluxM, fluxS )
    end
    
    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh, mesh2d );
    end
    
end

