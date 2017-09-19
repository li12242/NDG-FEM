classdef std_fv_tri < ndg_test.cell_test.std_fv_cell
    %STD_FV_TRI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Access=protected, Hidden=true)
        function [ NFV, EToV ] = loc_fv_info( obj )
        % set the vertex lists of the local finite volume element
            cell = obj.cell;
            NFV = cell.N^2;
            EToV = zeros(3, NFV);
            
            sk = 1;
            for row = 1:cell.N
                vert1 = cell.Fmask(row, 3);
                vert2 = cell.Fmask(row+1, 3);
                for kb = 1:row
                    EToV(:, sk) = [vert1, vert2, vert2+1]';
                    vert1 = vert1+1;
                    vert2 = vert2+1;
                    sk = sk+1;
                end

                vert1 = cell.Fmask(row, 3);
                vert2 = cell.Fmask(row+1, 3);
                for kt = 1:row-1
                    EToV(:, sk) = [vert1, vert2+1, vert1+1]';
                    vert1 = vert1+1;
                    vert2 = vert2+1;
                    sk = sk+1;
                end
            end
        end% func
        
        function [ Nedge, v1, v2 ] = loc_edge_info( obj )
            Nface = obj.cell.Nface;
            Nedge = obj.NFV*Nface;
            v1 = zeros(Nface, obj.NFV);
            v2 = zeros(Nface, obj.NFV);
            for k = 1:obj.NFV % loop over refined elements
                for f1 = 1:Nface
                    f2 = mod(f1,Nface)+1;
                    v1(f1, k) = obj.EToV(f1,k);
                    v2(f1, k) = obj.EToV(f2,k);
                end
            end
            v1 = v1(:); v2 = v2(:);
        end% func
        
        function [ vol ] = loc_cv_info(obj)
            vol = zeros(obj.NCV, 1);
            r = obj.cell.r;
            s = obj.cell.s;
            for k = 1:obj.NFV
                vert = obj.EToV(:, k);
                area = tri_area( r(vert), s(vert) );
                vol(vert) = vol(vert) + area/3;
            end
        end% func
        
        function [ P, R ] = loc_proj_recon_mat(obj)
            cell = obj.cell;
            Nface = cell.Nface;
            r = cell.r; s = cell.s;
            P = zeros( cell.Np, cell.Np );
            Pt = zeros( cell.Np, cell.Np);
            for k = 1:obj.NFV
                vert = obj.EToV(:, k);
                rc = mean( r(vert) );
                sc = mean( s(vert) );
                
                for f1 = 1:Nface
                    vind1 = f1;
                    vind2 = mod(vind1, cell.Nv)+1;
                    v1 = obj.EToV(vind1,k);
                    v2 = obj.EToV(vind2,k);
                    
                    r1 = r(v1); r2 = r(v2); rf = 0.5*(r1+r2);
                    s1 = s(v1); s2 = s(v2); sf = 0.5*(s1+s2);
                    
                    rq = cell.project_vert2node([rf, rc, r1]');
                    sq = cell.project_vert2node([sf, sc, s1]');
                    a1 = tri_area([r1, rf, rc]', [s1, sf, sc]')/2;
                    for i = 1:cell.Np
                        Pt(:,i) = cell.orthogonal_func(cell.N, i, rq, sq);
                    end
                    P(v1, :) = P(v1, :) + (a1.*cell.w')*Pt;
                    
                    rt2 = cell.project_vert2node([rf, r2, rc]');
                    st2 = cell.project_vert2node([sf, s2, sc]');
                    a2 = tri_area([r2, rf, rc]', [s2, sf, sc]')/2;
                    for i = 1:cell.Np
                        Pt(:,i) = cell.orthogonal_func(cell.N, i, rt2, st2);
                    end
                    P(v2, :) = P(v2, :) + (a2.*cell.w')*Pt;
                end
            end
            P = P/cell.V;
            P = bsxfun(@times, P, 1./obj.vol);
            R = inv(P);
        end% func
    end% methods
    
    methods
        function obj = std_fv_tri(cell)
            obj = obj@ndg_test.cell_test.std_fv_cell(cell);
        end% func
    end
    
end

function area = tri_area(vx, vy)
% calculate the area of triangle from given vertex coordinate
a = sqrt( (vx(1,:) - vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2 );
b = sqrt( (vx(2,:) - vx(3,:)).^2 + (vy(2,:)-vy(3,:)).^2 );
c = sqrt( (vx(3,:) - vx(1,:)).^2 + (vy(3,:)-vy(1,:)).^2 );
p = (a+b+c)./2;
area = sqrt(p.*(p-a).*(p-b).*(p-c));
end

