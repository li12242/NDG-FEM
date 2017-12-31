%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgEdge2d < NdgEdge
    
    methods
        function obj = NdgEdge2d( mesh1, mesh2, meshId1, meshId2 )
            obj = obj@NdgEdge( mesh1, mesh2, meshId1, meshId2 );
        end
        
%         function draw( obj )
%             vx = obj.umesh(1).vx( obj.FToV );
%             vy = obj.umesh(1).vy( obj.FToV );
%             plot(vx, vy, 'k.-');
%         end
    end
    
    methods(Hidden, Access=protected)
        function [ line1, line2 ] = setStdEdgeCell( obj, mesh1, mesh2 )
            N1 = mesh1.cell.N;
            N2 = mesh2.cell.N;
            line1 = StdLine( N1 );
            line2 = StdLine( N2 );
        end% func
        
%         function [ nxq, nyq, nzq, Jq ] = assembleQuadPointScale( obj, meshArray )
%             mesh = meshArray( obj.FToM(1) );
%             Nq = obj.bcell.Nq;
%             nxq = zeros(Nq, obj.M);
%             nyq = zeros(Nq, obj.M);
%             nzq = zeros(Nq, obj.M);
%             Jq = zeros(Nq, obj.M);
%             
%             for n = 1:obj.M
%                 vert = obj.FToV(:, n);
%                 vx = mesh.vx( vert );
%                 vy = mesh.vy( vert );
%                 nx =  ( vy(2) - vy(1) );
%                 ny = -( vx(2) - vx(1) );
%                 
%                 xc = mesh.xc( obj.FToE(1, n) );
%                 yc = mesh.yc( obj.FToE(1, n) );
%                 xfc = mean(vx);
%                 yfc = mean(vy);
%                 dx = xfc - xc(1);
%                 dy = yfc - yc(1);
%                 if( dx*nx + dy*ny ) < 0
%                     nx = -nx;
%                     ny = -ny;
%                 end
%                 
%                 J = sqrt(nx.*nx + ny.*ny);
%                 nxq(:, n) = nx./J;
%                 nyq(:, n) = ny./J;
%                 Jq(:, n) = J.*0.5;
%             end
%         end% func
    end
    
end

