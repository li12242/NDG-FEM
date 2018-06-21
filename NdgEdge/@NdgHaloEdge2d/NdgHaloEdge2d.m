%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgHaloEdge2d < NdgHaloEdge
    
    methods
        function obj = NdgHaloEdge2d( meshUnion, locMeshId, BCToV )
            obj = obj@NdgHaloEdge( meshUnion, locMeshId, BCToV );
        end
        
        [ fnode ] = proj_vert2node( obj, fvert );
        
        function draw(obj, varargin)
            ind = obj.FToN1 + ( obj.FToE - 1 ) * obj.mesh.cell.Np;
            % indK = repmat( obj.FToE(1, :), obj.Nfp, 1 );
            % indM = sub2ind( [obj.mesh.cell.Np, obj.mesh.K], obj.FToN1, indK );
            xv = obj.mesh.x( ind );
            yv = obj.mesh.y( ind );
            if nargin == 1
                plot( xv, yv, 'k.-', 'LineWidth', 2 );
            elseif nargin == 2
                plot3( xv, yv, varargin{1}, 'k.-', 'LineWidth', 2 );
            end
        end
    end
    
    methods( Access = protected )
        obj = assembleMassMatrix( obj );
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
        obj = assembleBoundaryConnection(obj, BCToV);
    end
    
end

