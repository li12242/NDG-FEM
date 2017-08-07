classdef tri_loc_fv2d < ndg_test.mesh_test.loc_fv
    %TRI_FV2D Finite volume information of the local refined element
    %   Detailed explanation goes here
    
    methods
        function obj = tri_loc_fv2d(mesh)
            obj = obj@ndg_test.mesh_test.loc_fv(mesh);
        end
        
        function draw(obj, varargin)
            hold on;
            patch('Vertices', [obj.mesh.cell.r(:), obj.mesh.cell.s(:)], ...
                'Faces', obj.EToV', ...
                'FaceColor', [0.8, 0.9, 1]);
            plot(obj.mesh.cell.r(:), obj.mesh.cell.s(:), 'ko', ...
                'MarkerFaceColor', 'k', ...
                'MarkerSize', 4)
        end
        
        function draw_control_vol(obj, n)
            obj.draw();
            r = obj.mesh.cell.r;
            s = obj.mesh.cell.s;
            rf = zeros(2, 1); 
            sf = zeros(2, 1);
            for k = 1:obj.Kloc
                vert = obj.EToV(:, k);
                fid = 1;
                if any(vert == n)
                    rc = mean( r(vert) ); sc = mean( s(vert) );
                    
                    vid = find(vert == n);
                    for t = 1:numel(vert)
                        if t == vid
                            continue;
                        end
                        rf(fid) = ( r(vert(t)) + r(n) )./2;
                        sf(fid) = ( s(vert(t)) + s(n) )./2;
                        fid = fid+1;
                    end
                    
                    patch('Vertices', ...
                        [r(n), rf(1), rc, rf(2); s(n), sf(1), sc, sf(2)]', ...
                        'Faces', [1,2,3,4], ...
                        'FaceColor', [0.9882, 0.9608, .8667], ...
                        'EdgeColor', 'k');
                    
                    plot([rf(1), rc, rf(2)], [sf(1), sc, sf(2)], 's-',...
                        'Color', 'k',...
                        'MarkerFaceColor', 'k', ...
                        'MarkerSize', 4)
                end
            end
            plot(r(n), s(n), 'ko-','Color', 'k',...
                'MarkerFaceColor', 'k', 'MarkerSize', 4)
        end
    end
    
    %% 私有方法
    methods(Access=protected)
        
        function [ Js ] = loc_fv_surface(obj, mesh)
            Js = zeros(size(mesh.Js)); % allocation
            
            Fmask = mesh.cell.Fmask;
            x = mesh.x; 
            y = mesh.y;
            for k = 1:mesh.K
                for n1 = 1:mesh.cell.Nfptotal-1
                    n2 = n1+1;
                    v1 = Fmask(n1);
                    v2 = Fmask(n2);
                    d2 = sqrt( (x(v1, k) - x(v2, k)).^2 ...
                        + (y(v1, k) - y(v2, k)).^2 )/2;
                    Js(n1, k) = Js(n1, k) + d2;
                    Js(n2, k) = Js(n2, k) + d2;
                end
            end
        end
        
        function [ P, R ] = project_matrix(obj, mesh)
            cell = mesh.cell;
            P = zeros(cell.Np, cell.Np);
            
            vol = zeros(cell.Np, 1);
            r = cell.r; 
            s = cell.s;
            xw = obj.TriGaussPoints(cell.N);
            w = xw(:, 3); w = w./sum(w); Nq = numel(w);
            Pw = zeros(Nq, cell.Np);
            for k = 1:obj.Kloc
                vert = obj.EToV(:, k);
                rc = mean( r( vert ) ); % cell centre coordinate
                sc = mean( s( vert ) );
                for f1 = 1:cell.Nface
                    v1 = f1;
                    v2 = mod(v1, cell.Nv)+1;
                    vert1 = obj.EToV(f1,k);
                    vert2 = obj.EToV(v2,k);
                    
                    r1 = r(vert1); s1 = s(vert1);
                    r2 = r(vert2); s2 = s(vert2);
                    rf = 0.5*(r1+r2);
                    sf = 0.5*(s1+s2);
                    
                    rt1 = obj.project_vert2quad(cell.N, [rf, rc, r1]');
                    st1 = obj.project_vert2quad(cell.N, [sf, sc, s1]');
                    a1 = obj.tri_area([r1, rf, rc]', [s1, sf, sc]');
                    vol(vert1) = vol(vert1) + a1;
                    for i = 1:cell.Np
                        Pw(:,i) = cell.orthogonal_func(cell.N, i, rt1, st1);
                    end
                    P(vert1, :) = P(vert1, :) + (a1.*w')*Pw;
                    
                    rt2 = obj.project_vert2quad(cell.N, [rf, r2, rc]');
                    st2 = obj.project_vert2quad(cell.N, [sf, s2, sc]');
                    a2 = obj.tri_area([r2, rf, rc]', [s2, sf, sc]');
                    vol(vert2) = vol(vert2) + a2;
                    for i = 1:cell.Np
                        Pw(:,i) = cell.orthogonal_func(cell.N, i, rt2, st2);
                    end
                    P(vert2, :) = P(vert2, :) + (a2.*w')*Pw;
                end
            end% for
            P = bsxfun(@times, P, 1./vol);
            P = P/cell.V;
            R = inv(P);
        end
        
        function [v1, v2, nx, ny, nz, ds] = loc_edge_info(obj, mesh)
            % find vertex pairs in each sub-cell
            K = mesh.K;
            Np = mesh.cell.Np;
            Nface = mesh.cell.Nface;
            v1 = zeros(obj.Nedge, K);
            v2 = zeros(obj.Nedge, K);
            nx = zeros(obj.Nedge, K);
            ny = zeros(obj.Nedge, K);
            nz = zeros(obj.Nedge, K);
            ds = zeros(obj.Nedge, K);
            
            xc = mesh.cell_mean(mesh.x);
            yc = mesh.cell_mean(mesh.y);
            for e = 1:K
                sk = 1;
                Npoint = (e-1)*Np;
                for k = 1:obj.Kloc % loop over refined elements
                    for f1 = 1:Nface
                        f2 = mod(f1,Nface)+1;
                        v1(sk, e) = obj.EToV(f1,k) + Npoint;
                        v2(sk, e) = obj.EToV(f2,k) + Npoint;
                        sk = sk + 1;
                    end
                end
                [ nx(:, e), ny(:, e), ds(:, e) ] = obj.norm_vector ...
                    ( v1(:, e), v2(:, e), mesh.x, mesh.y, xc(e), yc(e) );
            end
        end% func
                
        function [Nedge, Kloc, EToV] = loc_connect(obj, cell)
            % 单元内
            Kloc = cell.N^2;
            Nedge = Kloc*cell.Nface;
            EToV = zeros(3, Kloc);
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
        
        function [ vol ] = loc_fv_info(obj, mesh)
            K = mesh.K;
            vol = zeros(mesh.cell.Np, K);
            for e = 1:K
                for k = 1:obj.Kloc
                    v = obj.EToV(:, k);
                    x = mesh.x(v, e);
                    y = mesh.y(v, e);
                    a = obj.tri_area(x, y);
                    vol(v, e) = vol(v, e) + a/3;
                end
            end
        end% func
    end
    
    %% private function
    methods(Access=private)
        function node_val = project_vert2quad(obj, N, vert_val)
            xw = obj.TriGaussPoints(N);
            r = xw(:, 1); % standard tri [0,1]x[0,1]
            s = xw(:, 2);
            node_val = ( (1-r-s)*vert_val(1, :) ...
                + r*vert_val(2, :) + s*vert_val(3, :) );
        end% func
    end
  
    methods(Access=private, Static)
        
        function [nx, ny, ds] = norm_vector(v1, v2, x, y, xc, yc)
            Nedge = numel(v1);
            nx = zeros(Nedge, 1);
            ny = zeros(Nedge, 1);
            ds = zeros(Nedge, 1);
            for f = 1:Nedge
                vert1 = v1(f);
                vert2 = v2(f);
                xv1 = x(vert1); yv1 = y(vert1); % vertex 
                xv2 = x(vert2); yv2 = y(vert2);
                xfc = 0.5*(xv1 + xv2); % middle points coordinate
                yfc = 0.5*(yv1 + yv2);
                
                nxt =  (yfc - yc);
                nyt = -(xfc - xc);
                ds(f) = sqrt( nxt*nxt + nyt*nyt );
                % check the direction is from v1 to v2
                dx = xv2 - xv1;
                dy = yv2 - yv1;
                if ( dx*nxt + dy*nyt ) < 0
                    nx(f) = -nxt/ds(f);
                    ny(f) = -nyt/ds(f);
                else
                    nx(f) = nxt/ds(f);
                    ny(f) = nyt/ds(f);
                end
            end
        end
                
        function area = tri_area(x, y)
            % calculate the area of triangle from given points coordinate
            a = sqrt( (x(1,:) - x(2,:)).^2 + (y(1,:)-y(2,:)).^2 );
            b = sqrt( (x(2,:) - x(3,:)).^2 + (y(2,:)-y(3,:)).^2 );
            c = sqrt( (x(3,:) - x(1,:)).^2 + (y(3,:)-y(1,:)).^2 );
            p = (a+b+c)./2;
            area = sqrt(p.*(p-a).*(p-b).*(p-c));
        end
                
        function xw = TriGaussPoints(n)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function TriGaussPoints provides the Gaussian points and weights %
        % for the Gaussian quadrature of order n for the standard triangles. %
        % http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF %
        % Input: n - the order of the Gaussian quadrature (n<=12) %
        % %
        % Output: xw - a n by 3 matrix: %
        % 1st column gives the x-coordinates of points %
        % 2nd column gives the y-coordinates of points %
        % 3rd column gives the weights %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (n == 1)
            xw=[0.33333333333333 0.33333333333333 1.00000000000000];
            elseif (n == 2)
            xw=[0.16666666666667 0.16666666666667 0.33333333333333
                0.16666666666667 0.66666666666667 0.33333333333333
                0.66666666666667 0.16666666666667 0.33333333333333];
            elseif (n == 3)
            xw=[0.33333333333333 0.33333333333333 -0.56250000000000
                0.20000000000000 0.20000000000000 0.52083333333333
                0.20000000000000 0.60000000000000 0.52083333333333
                0.60000000000000 0.20000000000000 0.52083333333333];
            elseif (n == 4)
            xw=[0.44594849091597 0.44594849091597 0.22338158967801
                0.44594849091597 0.10810301816807 0.22338158967801
                0.10810301816807 0.44594849091597 0.22338158967801
                0.09157621350977 0.09157621350977 0.10995174365532
                0.09157621350977 0.81684757298046 0.10995174365532
                0.81684757298046 0.09157621350977 0.10995174365532];
            elseif (n == 5)
            xw=[0.33333333333333 0.33333333333333 0.22500000000000
                0.47014206410511 0.47014206410511 0.13239415278851
                0.47014206410511 0.05971587178977 0.13239415278851
                0.05971587178977 0.47014206410511 0.13239415278851
                0.10128650732346 0.10128650732346 0.12593918054483
                0.10128650732346 0.79742698535309 0.12593918054483
                0.79742698535309 0.10128650732346 0.12593918054483];
            elseif (n == 6)
            xw=[0.24928674517091 0.24928674517091 0.11678627572638
                0.24928674517091 0.50142650965818 0.11678627572638
                0.50142650965818 0.24928674517091 0.11678627572638
                0.06308901449150 0.06308901449150 0.05084490637021
                0.06308901449150 0.87382197101700 0.05084490637021
                0.87382197101700 0.06308901449150 0.05084490637021
                0.31035245103378 0.63650249912140 0.08285107561837
                0.63650249912140 0.05314504984482 0.08285107561837
                0.05314504984482 0.31035245103378 0.08285107561837
                0.63650249912140 0.31035245103378 0.08285107561837
                0.31035245103378 0.05314504984482 0.08285107561837
                0.05314504984482 0.63650249912140 0.08285107561837];
            elseif (n == 7)
            xw=[0.33333333333333 0.33333333333333 -0.14957004446768
                0.26034596607904 0.26034596607904 0.17561525743321
                0.26034596607904 0.47930806784192 0.17561525743321
                0.47930806784192 0.26034596607904 0.17561525743321
                0.06513010290222 0.06513010290222 0.05334723560884
                0.06513010290222 0.86973979419557 0.05334723560884
                0.86973979419557 0.06513010290222 0.05334723560884
                0.31286549600487 0.63844418856981 0.07711376089026
                0.63844418856981 0.04869031542532 0.07711376089026
                0.04869031542532 0.31286549600487 0.07711376089026
                0.63844418856981 0.31286549600487 0.07711376089026
                0.31286549600487 0.04869031542532 0.07711376089026
                0.04869031542532 0.63844418856981 0.07711376089026];
            elseif (n == 8)
            xw=[0.33333333333333 0.33333333333333 0.14431560767779
                0.45929258829272 0.45929258829272 0.09509163426728
                0.45929258829272 0.08141482341455 0.09509163426728
                0.08141482341455 0.45929258829272 0.09509163426728
                0.17056930775176 0.17056930775176 0.10321737053472
                0.17056930775176 0.65886138449648 0.10321737053472
                0.65886138449648 0.17056930775176 0.10321737053472
                0.05054722831703 0.05054722831703 0.03245849762320
                0.05054722831703 0.89890554336594 0.03245849762320
                0.89890554336594 0.05054722831703 0.03245849762320
                0.26311282963464 0.72849239295540 0.02723031417443
                0.72849239295540 0.00839477740996 0.02723031417443
                0.00839477740996 0.26311282963464 0.02723031417443
                0.72849239295540 0.26311282963464 0.02723031417443
                0.26311282963464 0.00839477740996 0.02723031417443
                0.00839477740996 0.72849239295540 0.02723031417443];
            else
                error('Bad input n');
            end
            return;
        end% func
    end
    
end

