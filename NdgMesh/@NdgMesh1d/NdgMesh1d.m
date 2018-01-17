classdef NdgMesh1d < NdgMesh
    %LINE_MESH Summary of this class goes here
    %   Detailed explanation goes here

    properties(Constant)
        type = NdgMeshType.OneDim
    end
    
    methods(Hidden, Access = protected)
        function [ edge ] = makeConnectNdgEdge( obj, mesh1, mid0, mid1 )
            edge = NdgEdge1d( obj, mesh1, mid0, mid1 );
        end
        
        function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J ] = assembleJacobiFactor( obj )
            xr = obj.cell.Dr*obj.x;
            J = xr; rx = 1./J;

            ry = zeros(size(rx));
            rz = zeros(size(rx)); 
            sx = zeros(size(rx));
            sy = zeros(size(rx));
            sz = zeros(size(rx)); 
            tx = zeros(size(rx));
            ty = zeros(size(rx));
            tz = zeros(size(rx)); 
        end
        
        function [nx, ny, nz, Js] = assembleFacialJaobiFactor( obj )
            % Define outward normals
            xb = obj.x(obj.cell.Fmask, :);
            xc = obj.GetMeshAverageValue( obj.x );
            nx = sign( xb - xc );
            
            Js = ones(size(nx));
            ny = zeros(obj.cell.Nface, obj.K);
            nz = zeros(obj.cell.Nface, obj.K);
        end
        
        function faceId = assembleGlobalFaceIndex( obj )
            faceId = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                faceId(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
            end
        end% func
    end% methods
    
    methods
        obj = refine(obj, refine_level);
        
        function obj = NdgMesh1d( cell, Nv, vx, K, EToV, EToR, BCToV )
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros            
            obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV);
        end% func
        
        
        function draw( obj, zvar )
            % check the figure is created and not removed.
            isFigureExit = ~isempty( obj.figureHandle ) ...
                && isvalid( obj.figureHandle ) ;
            
            if isFigureExit
                %updateFigure(obj.figureHandle, obj.x(:), zvar(:));
                set( obj.figureHandle, 'YData', zvar(:) );
            else
                Np = obj.cell.Np; 
                K = obj.K;
                list = 1:Np:(K*Np);
                g = graph();
                for n = 1:Np-1
                    g = addedge(g, list+n-1, list+n);
                end
                
                obj.figureHandle = plot(g, ...
                    'XData', obj.x(:), 'YData', zvar(:), ...
                    'LineWidth', 1, ...
                    'Marker', 'o', ...
                    'NodeColor','k', ...
                    'EdgeColor', 'k', ...
                    'MarkerSize', 2, ...
                    'NodeLabel', {});
                box on;
                grid on;
            end

        end
    end
    
end

