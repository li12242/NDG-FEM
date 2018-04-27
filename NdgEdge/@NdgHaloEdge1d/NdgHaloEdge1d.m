classdef NdgHaloEdge1d < NdgHaloEdge
    
    methods
        function obj = NdgHaloEdge1d( mesh1, mesh2, meshId1, meshId2 )
            obj = obj@NdgHaloEdge( mesh1, mesh2, meshId1, meshId2 );
        end
        
%         function draw( obj )
%             plot()
%         end
    end
    
    methods(Hidden, Access=protected)
        
        function [ point1, point2 ] = setStdEdgeCell( obj, mesh1, mesh2 )
            N1 = mesh1.cell.N;
            N2 = mesh2.cell.N;
            point1 = StdPoint( N1 );
            point2 = StdPoint( N2 );
        end
        
        %> Calculate the outward normal vector and the determination of
        %> Jacobian matrix at quadrature points on each edge.
%         function [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale( obj )
%             mesh = obj.umesh(1);
%             Nq = obj.bcell.Nq;
%             nxq = ones(Nq, obj.M);
%             nyq = zeros(Nq, obj.M);
%             nzq = zeros(Nq, obj.M);
%             Jq = ones(Nq, obj.M);
%             
%             xc = mean( mesh.x( :, obj.FToE(1,:) ) );
%             % Define outward normals
%             for n = 1:obj.M
%                 xb = mesh.vx( obj.FToV(:, n) );
%                 if ( xb - xc(n) ) < 1e-10
%                     nxq(n) = -1;
%                 end
%             end
%         end% func
    end
    
end

