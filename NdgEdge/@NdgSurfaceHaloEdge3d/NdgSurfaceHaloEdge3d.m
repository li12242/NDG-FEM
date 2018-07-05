classdef NdgSurfaceHaloEdge3d < NdgBottomHaloEdge3d
    %NDGSURFACEHALOEDGE3D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = NdgSurfaceHaloEdge3d( meshUnion3d, meshId )
            obj = obj@NdgBottomHaloEdge3d( meshUnion3d, meshId );
        end
    end
    
    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh, mesh2d );
    end
    
end

