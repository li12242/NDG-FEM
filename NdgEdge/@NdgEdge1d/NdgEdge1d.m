classdef NdgEdge1d < NdgEdge
    
    methods
        function obj = NdgEdge1d(umesh, meshId, BCToV)
            obj = obj@NdgEdge(umesh, meshId, BCToV);
        end
        
        function draw(obj)
            plot()
        end
    end
    
    methods(Hidden, Access=protected)
        
        function [bcell] = setStdEdgeCell(obj)
            bcell = StdPoint(0);
        end
        
        %> Calculate the outward normal vector and the determination of
        %> Jacobian matrix at quadrature points on each edge.
        function [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale( obj )
            mesh = obj.umesh(1);
            Nq = obj.bcell.Nq;
            nxq = ones(Nq, obj.M);
            nyq = zeros(Nq, obj.M);
            nzq = zeros(Nq, obj.M);
            Jq = ones(Nq, obj.M);
            
            xc = mean( mesh.x( :, obj.FToE(1,:) ) );
            % Define outward normals
            for n = 1:obj.M
                xb = mesh.vx( obj.FToV(:, n) );
                if ( xb - xc(n) ) < 1e-10
                    nxq(n) = -1;
                end
            end
        end% func
    end
    
end

