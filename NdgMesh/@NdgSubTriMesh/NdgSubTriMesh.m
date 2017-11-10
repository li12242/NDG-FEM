classdef NdgSubTriMesh < NdgMesh2d
    %NDGSUBTRIMESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        locNx, locNy, locNz
        locEdgeJs
        LAVToCV
    end
    
    methods
        
        function obj = NdgSubTriMesh(cell, Nv, vx, vy, K, EToV, EToR, BCToV)
            obj = obj@NdgMesh2d(cell, Nv, vx, vy, K, EToV, EToR, BCToV);
            [ obj.locNx, obj.locNy, obj.locNz, obj.locEdgeJs ] = assembleLocalEdgeInfo( obj );
            [ obj.LAVToCV ] = assembleSubCellLAV( obj );
            [ obj.Js ] = assembleCellEdgeInfo( obj );
        end
    end
    
    methods( Hidden, Access = protected )
        function [ nx, ny, nz, locJs ] = assembleLocalEdgeInfo( obj )
            xr = obj.cell.Dr * obj.x;
            yr = obj.cell.Dr * obj.y;
            xs = obj.cell.Ds * obj.x;
            ys = obj.cell.Ds * obj.y;
            
            fxr = xr(1,:); fyr = yr(1,:);
            fxs = xs(1,:); fys = ys(1,:);
            
            nx = bsxfun(@times,  fys, obj.cell.nr) ...
                + bsxfun(@times, -fyr, obj.cell.ns);
            ny = bsxfun(@times, -fxs, obj.cell.nr) ...
                + bsxfun(@times,  fxr, obj.cell.ns);
            J  = sqrt(nx.^2 + ny.^2);
            nx = nx ./ J;
            ny = ny ./ J;
            nz = zeros( obj.cell.Nedge, obj.K );
            
            locJs = zeros( obj.cell.Nedge, obj.K );
            Nface = obj.cell.Nface;
            for k = 1:obj.K
                for n = 1:obj.cell.NFV
                    vind = obj.cell.EToV(:, n);
                    xc = mean( obj.x(vind, k) );
                    yc = mean( obj.y(vind, k) );
                    for f = 1:Nface
                        edgeId = (n-1)*Nface + f; % edge index
                        n1 = obj.cell.n1(edgeId);
                        n2 = obj.cell.n2(edgeId);
                        x1 = obj.x(n1, k);
                        y1 = obj.y(n1, k);
                        x2 = obj.x(n2, k);
                        y2 = obj.y(n2, k);
                        xfc = 0.5*(x1 + x2);
                        yfc = 0.5*(y1 + y2);
                        
                        nxt =  ( yfc - yc );
                        nyt = -( xfc - xc );
                        locJs(edgeId, k) = sqrt( nxt*nxt + nyt*nyt );
                    end
                end
            end
        end% func
        
        function [ LAVToCV ] = assembleSubCellLAV( obj )
            LAVToCV = bsxfun(@times, obj.cell.LAVToCV, obj.LAV)./2;
        end
        
        function [ Js ] = assembleCellEdgeInfo( obj )
            Js = zeros( size(obj.Js) );
            Fmask = obj.cell.Fmask;
            x = obj.x;
            y = obj.y;
            for k = 1:obj.K
                for n1 = 1:(obj.cell.TNfp-1)
                    n2 = n1 + 1;
                    v1 = Fmask(n1);
                    v2 = Fmask(n2);
                    d2 = sqrt( (x(v1, k) - x(v2, k)).^2 ...
                        + (y(v1, k) - y(v2, k)).^2 )/2;
                    Js(n1, k) = Js(n1, k) + d2;
                    Js(n2, k) = Js(n2, k) + d2;
                end
            end
            %Js = Js ./ obj.LAVToCV( obj.cell.Fmask(:), : );
        end% func
    end
    
end
