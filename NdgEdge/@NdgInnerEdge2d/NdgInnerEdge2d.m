classdef NdgInnerEdge2d < NdgInnerEdge
    
    methods
        function obj = NdgInnerEdge2d( mesh, meshId )
            obj = obj@NdgInnerEdge( mesh, meshId );
        end
        
        function draw(obj, varargin)
            xfM = zeros( obj.Nfp, obj.Ne );
            yfM = zeros( obj.Nfp, obj.Ne );
            mesh = obj.mesh;
            for i = 1:obj.Ne
                k1 = obj.FToE(1, i);                
                fp1 = obj.FToN1(:, i);
                xfM(:, i) = mesh.x( fp1, k1 );
                yfM(:, i) = mesh.y( fp1, k1 );
            end
            if (nargin == 2)
                plot3( xfM, yfM, varargin{1}, 'k.-', 'LineWidth', 2 );
            elseif (nargin == 1)
                plot( xfM, yfM, 'k.-', 'LineWidth', 2 );
            end
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

