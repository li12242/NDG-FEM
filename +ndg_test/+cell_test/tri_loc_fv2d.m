classdef tri_loc_fv2d
    %TRI_FV2D Finite volume information of the local refined element
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        Nedge % number of edges in each refined element
        Kloc % number of sub-cells in each refined element
        v1, v2 % vertex list of each edge
        nx, ny % normal vector pointing from v1 to v2
        ds % length of boundary edge    
        P, R % project and reconstruct matrix from lagrange basis 
             % coefficients to finite volume values, satisfying R = P^-1
    end
    
    methods
        function obj = tri_loc_fv2d(mesh)
            cell = mesh.cell;
            % connection in refined triangle
            [ obj.Kloc, EToV ] = obj.tri_local_connect(cell); 
            obj.Nedge = obj.Kloc*mesh.cell.Nface;
            
            [ obj.v1, obj.v2 ] = loc_fv_info(obj, mesh.K, ...
                cell.Nface, cell.Np, EToV);
            [ obj.nx, obj.ny, obj.ds ] = loc_fv_scale(obj, ...
                cell.Nface, mesh.K, EToV, mesh.x, mesh.y);
            [ obj.P, obj.R ] = project_matrix(obj, EToV, cell);
        end
        
        function val = project_node2fvm(obj, f_Q)
            val = obj.P*f_Q;
        end
        
        function f_Q = project_fvm2node(obj, val)
            f_Q = obj.R*val;
        end
    end
    
    methods(Access=protected)
        
        function [P, R, vol] = project_matrix(obj, EToV, cell)
            P = zeros(cell.Np, cell.Np);
            R = zeros(cell.Np, cell.Np);
            vol = zeros(cell.Np, 1);
            r = cell.r; 
            s = cell.s;
            
            xw = obj.TriGaussPoints(cell.N);
            w = xw(:, 3); Nq = numel(w);
            Pw = zeros(Nq, cell.Np);
            for k = 1:obj.Kloc
                vert = EToV(:, k);
                rc = mean( r( vert ) ); % cell centre coordinate
                sc = mean( s( vert ) );
                for f1 = 1:cell.Nface
                    f2 = mod(f1,cell.Nface)+1;
                    vert1 = EToV(f1,k);
                    vert2 = EToV(f2,k);
                    
                    r1 = r(vert1); s1 = s(vert1);
                    r2 = r(vert2); s2 = s(vert2);
                    rf = 0.5*(r1+r2);
                    sf = 0.5*(s1+s2);
                    
                    rt1 = obj.project_vert2quad(cell.N, [rf, rc, r1]');
                    st1 = obj.project_vert2quad(cell.N, [sf, sc, s1]');
                    a1 = obj.tri_area([r1, s1], [rf, sf], [rc, sc]);
                    vol(vert1) = vol(vert1) + a1;
                    for i = 1:cell.Np
                        Pw(:,i) = cell.orthogonal_func(cell.N, i, rt1, st1);
                    end
                    P(vert1, :) = P(vert1, :) + (a1.*w')*Pw;
                    
                    rt2 = obj.project_vert2quad(cell.N, [rf, r2, rc]');
                    st2 = obj.project_vert2quad(cell.N, [sf, s2, sc]');
                    a2 = obj.tri_area([r1, s1], [rf, sf], [rc, sc]);
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
        
        function [v1, v2] = loc_fv_info(obj, K, Nface, Np, EToV)
            % find vertex pairs in each sub-cell
            v1 = zeros(obj.Nedge, K);
            v2 = zeros(obj.Nedge, K);
            for e = 1:K
                sk = 1;
                Npoint = (e-1)*Np;
                for k = 1:obj.Kloc % loop over refined elements
                    for f1 = 1:Nface
                        f2 = mod(f1,Nface)+1;
                        v1(sk, e) = EToV(f1,k) + Npoint;
                        v2(sk, e) = EToV(f2,k) + Npoint;
                        sk = sk + 1;
                    end
                end% for
            end
        end
        
        function [nx, ny, ds] = loc_fv_scale(obj, Nface, K, EToV, x, y)
            nx = zeros(obj.Nedge, K);
            ny = zeros(obj.Nedge, K);
            ds = zeros(obj.Nedge, K);
            %xf = zeros(obj.Nedge, K);
            %yf = zeros(obj.Nedge, K);

            for e = 1:K % loop over all elements
                sk = 1;
                for k = 1:obj.Kloc
                    vert = EToV(:, k);
                    x1 = mean( x( vert, e ) ); % cell centre coordinate
                    y1 = mean( y( vert, e ) );
                    for f1 = 1:Nface
                        vert1 = obj.v1( sk, e );
                        vert2 = obj.v2( sk, e );
                        xv1 = x(vert1); yv1 = y(vert1); % vertex 
                        xv2 = x(vert2); yv2 = y(vert2);

                        x2 = 0.5*(xv1 + xv2); % middle points coordinate
                        y2 = 0.5*(yv1 + yv2);
                        %xf(sk, e) = x2;
                        %yf(sk, e) = y2;

                        nxt =  (y2 - y1);
                        nyt = -(x2 - x1);
                        ds(sk, e) = sqrt( nxt*nxt + nyt*nyt );
                        % check the direction is from v1 to v2
                        dx = xv2 - xv1;
                        dy = yv2 - yv1;
                        if ( dx*nxt + dy*nyt ) < 0
                            nx(sk, e) = -nxt/ds(sk, e);
                            ny(sk, e) = -nyt/ds(sk, e);
                        else
                            nx(sk, e) = nxt/ds(sk, e);
                            ny(sk, e) = nyt/ds(sk, e);
                        end
                        sk = sk + 1;
                    end                
                end
            end
            %quiver(xf, yf, nx.*ds, ny.*ds);
        end
        
        function node_val = project_vert2quad(obj, N, vert_val)
            xw = obj.TriGaussPoints(N);
            r = xw(:, 1); 
            s = xw(:, 2);
            node_val = ( (1-r-s)*vert_val(1, :) ...
                + r*vert_val(2, :) + s*vert_val(3, :) );
        end% 
    end
    
    methods(Access=private, Static)
        function area = tri_area(p1, p2, p3)
            a = sqrt( (p1(1)-p2(1))^2+(p1(2)-p2(2))^2 );
            b = sqrt( (p3(1)-p2(1))^2+(p3(2)-p2(2))^2 );
            c = sqrt( (p1(1)-p3(1))^2+(p1(2)-p3(2))^2 );
            p = (a+b+c)/2;
            area = sqrt(p*(p-a)*(p-b)*(p-c));
        end
        
        function [Kloc, EToV] = tri_local_connect(cell)
            Kloc = cell.N^2;
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

