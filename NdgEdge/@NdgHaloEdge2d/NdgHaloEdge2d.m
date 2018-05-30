%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgHaloEdge2d < NdgHaloEdge
    
    methods
        function obj = NdgHaloEdge2d( meshUnion, locMeshId )
            obj = obj@NdgHaloEdge( meshUnion, locMeshId );
        end
        
        function draw(obj, varargin)
            indK = repmat( obj.FToE(1, :), obj.Nfp, 1 );
            indM = sub2ind( [obj.mesh.cell.Np, obj.mesh.K], obj.FToN1, indK );
            xv = obj.mesh.x( indM );
            yv = obj.mesh.y( indM );
            if nargin == 1
                plot( xv, yv, 'k.-', 'LineWidth', 2 );
            elseif nargin == 2
                plot3( xv, yv, varargin{1}, 'k.-', 'LineWidth', 2 );
            end
        end
    end
    
    methods(Hidden = true, Access = protected)
        function [ eCell ] = setEdgeReferCell( obj, mesh )
            eCell = StdLine( mesh.cell.N );
        end
    end
    
    methods( Access = protected )
        [ Nedge, FToE, FToF, FToV, FToM, ftype ] = assembleEdgeConnect( obj, mesh, meshId );
        [ FToN1, FToN2, nx, ny, nz, Js ] = assembleNodeProject( obj, meshUnion )
    end
    
end

