classdef mesh_loc_tri < ndg_test.mesh_test.mesh_loc_fv
    %MESH_LOC_TRI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = mesh_loc_tri(subcell, mesh)
            obj = obj@ndg_test.mesh_test.mesh_loc_fv(subcell, mesh);
        end
    end
    
    methods(Hidden=true, Access=protected)
        function [ nx, ny, nz, ds ] = inner_edge_info(obj)
            Nface = 3;
            nx = zeros(obj.subcell.Nedge, obj.mesh.K);
            ny = zeros(obj.subcell.Nedge, obj.mesh.K);
            nz = zeros(obj.subcell.Nedge, obj.mesh.K);
            ds = zeros(obj.subcell.Nedge, obj.mesh.K);
             
            mesh = obj.mesh;
            for k = 1:mesh.K
                for e = 1:obj.subcell.NFV
                    vind = obj.subcell.EToV(:, e);
                    xc = mean( mesh.x(vind, k) );
                    yc = mean( mesh.y(vind, k) );
                    for f = 1:Nface
                        n = (e-1)*Nface + f; % edge index
                        v1 = obj.subcell.v1(n);
                        v2 = obj.subcell.v2(n);
                        x1 = mesh.x(v1, k); y1 = mesh.y(v1, k);
                        x2 = mesh.x(v2, k); y2 = mesh.y(v2, k);
                        xfc = 0.5*(x1 + x2);
                        yfc = 0.5*(y1 + y2);

                        nxt =  ( yfc - yc );
                        nyt = -( xfc - xc );
                        ds(n, k) = sqrt( nxt*nxt + nyt*nyt );
                        % check the direction is from v1 to v2
                        dx = x2 - x1;
                        dy = y2 - y1;
                        if ( dx*nxt + dy*nyt ) < 0
                            nx(n, k) = -nxt/ds(n, k);
                            ny(n, k) = -nyt/ds(n, k);
                        else
                            nx(n, k) = nxt/ds(n, k);
                            ny(n, k) = nyt/ds(n, k);
                        end
                    end
                end
            end
        end% func
        
        function [ vol ] = loc_cv_vol(obj)
            mesh = obj.mesh;
            vol = zeros(obj.subcell.NCV, mesh.K);
            for k = 1:mesh.K
                for n = 1:obj.subcell.NFV
                    vind = obj.subcell.EToV(:, n);
                    vx = mesh.x(vind, k);
                    vy = mesh.y(vind, k);
                    area = tri_area(vx, vy);
                    vol(vind, k) = vol(vind, k) + area/3;
                end
            end
        end
        
        function [ Js ] = cv_surf_info(obj)
            mesh = obj.mesh;
            Js = zeros( size(mesh.Js) ); % allocation
            
            Fmask = mesh.cell.Fmask;
            x = mesh.x; 
            y = mesh.y;
            for k = 1:mesh.K
                for n1 = 1:(mesh.cell.Nfptotal-1)
                    n2 = n1 + 1;
                    v1 = Fmask(n1);
                    v2 = Fmask(n2);
                    d2 = sqrt( (x(v1, k) - x(v2, k)).^2 ...
                        + (y(v1, k) - y(v2, k)).^2 )/2;
                    Js(n1, k) = Js(n1, k) + d2;
                    Js(n2, k) = Js(n2, k) + d2;
                end
            end
        end
    end% methods
end

function area = tri_area(vx, vy)
% calculate the area of triangle from given vertex coordinate
a = sqrt( (vx(1,:) - vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2 );
b = sqrt( (vx(2,:) - vx(3,:)).^2 + (vy(2,:)-vy(3,:)).^2 );
c = sqrt( (vx(3,:) - vx(1,:)).^2 + (vy(3,:)-vy(1,:)).^2 );
p = (a+b+c)./2;
area = sqrt(p.*(p-a).*(p-b).*(p-c));
end

