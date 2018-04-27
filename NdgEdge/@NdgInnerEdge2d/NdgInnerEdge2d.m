classdef NdgInnerEdge2d < NdgInnerEdge
    
    methods
        function obj = NdgInnerEdge2d( mesh, meshId )
            obj = obj@NdgInnerEdge( mesh, meshId );
        end        
    end
    
    methods( Static, Access = protected )
        %> set reference standard cell
        function [ bcell ] = setEdgeReferCell( mesh )
            bcell = StdLine( mesh.cell.N );
        end
    end
    
end

