classdef NdgInnerEdge2d < NdgInnerEdge
    
    methods
        function obj = NdgInnerEdge2d( mesh, meshId )
            obj = obj@NdgInnerEdge( mesh, meshId );
        end
        
        function draw(obj)
            xv = obj.mesh.vx( obj.FToV );
            yv = obj.mesh.vy( obj.FToV );
            plot( xv, yv, 'k.-', 'LineWidth', 2 );
        end
    end
    
    methods( Static, Access = protected )
        %> set reference standard cell
        function [ eCell ] = setEdgeReferCell( mesh )
            eCell = StdLine( mesh.cell.N );
        end
    end
    
    methods( Access = protected )
        %> connect edge to elements
        [ Nedge, FToE, FToF, FToV ] = assembleEdgeConnect( obj, mesh )
        [ FToN1, FToN2, nx, ny, nz, Js ] = assembleNodeProject( obj, mesh )
    end
    
end

