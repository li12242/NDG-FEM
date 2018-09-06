classdef NdgInnerEdge1d < NdgInnerEdge
    %NDGINNEREDGE1D Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = NdgInnerEdge1d( meshUnion, meshId )
            obj = obj@NdgInnerEdge( meshUnion, meshId );
        end
    end
    
    methods( Access = protected )
        obj = assembleMassMatrix( obj );
        %> connect edge to elements
        obj = assembleEdgeConnect( obj, mesh )
        obj = assembleNodeProject( obj, mesh )
    end
end

