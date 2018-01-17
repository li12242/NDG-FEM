classdef StdSubCellTest < matlab.unittest.TestCase
    %STDSUBCELLTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(MethodSetupParameter)
        %> test cell types
        type = {...
            NdgCellType.Line, ...
            NdgCellType.Tri, ...
%             NdgCellType.Quad
            }
        %> test cell orders
        order = {1,2,3,4}
    end
    
    properties(Constant)
        %> tolerance
        tol = 1e-9;
    end
    
    properties
        %> cell object
        cell
    end
    
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_std_cell(test, type, order)
            switch type
                case NdgCellType.Line
                    test.cell = StdSubLine( order );
                case NdgCellType.Tri
                    test.cell = StdSubTri( order );
            end
        end% func
    end
    
    methods(Test, ParameterCombination = 'sequential')
        function TestEdgeNormalVector( test )
%             stdcell = test.cell;
%             [nx, ny, nz, flen] = stdcell.matEvaluateInnerControlEdgeInfo...
%                 ( stdcell.r, stdcell.s, stdcell.t );
%             % evaluate the coordinate at the edge
%             E = 0.5.* ( stdcell.r( stdcell.n1 ) + stdcell.r( stdcell.n2 ) );
%             G = 0.5.* ( stdcell.s( stdcell.n1 ) + stdcell.s( stdcell.n2 ) );
%             H = 0.5.* ( stdcell.t( stdcell.n1 ) + stdcell.t( stdcell.n2 ) );
%             flux = flen.*(E.*nx + G.*ny + H.*nz);
%             rhs = (stdcell.VLIFT * flux);
%             
%             E = stdcell.r( stdcell.Fmask(:) );
%             G = stdcell.r( stdcell.Fmask(:) );
%             H = stdcell.r( stdcell.Fmask(:) );
%             [ nx, ny, nz, Js ] = stdcell.assembleNormalVector...
%                 ( stdcell.r, stdcell.s, stdcell.t );
%             flux = Js.*(E.*nx + G.*ny + H.*nz);
%             rhs = rhs + stdcell.SLIFT * flux;
%             
%             rhs = rhs ./ stdcell.LAVToCV;
%             rhs_ext = ones( size(rhs) );
%             test.verifyEqual(rhs, rhs_ext, 'AbsTol', test.tol);
        end
        
        function TestProjectMatrix( test )
            stdcell = test.cell;
            switch stdcell.type
                case NdgCellType.Line
                    TestLineProjectMatrix( test );
                case NdgCellType.Tri
                    TestTriProjectMatrix( test );
                otherwise
            end
        end
    end% methods
    
    methods
        function TestLineProjectMatrix( test )
            stdcell = test.cell;
            r = stdcell.r;
            FunHandle = @( r,n ) ( r ).^(n);
            for n = 0:stdcell.N
                f_ext = zeros( stdcell.Np, 1 );
                for i = 1:stdcell.NFV
                    vert = stdcell.EToV(:, i);
                    rv = r( vert );
                    rc = mean( rv ); 
                    area = abs( diff( rv ) );
                    
                    rv1 = [rv(1), rc]';
                    rv2 = [rc, rv(2)]';
                    rq1 = stdcell.project_vert2quad(rv1);
                    rq2 = stdcell.project_vert2quad(rv2);
                    fval = area./2 .* ( stdcell.wq' ...
                        * FunHandle( rq1, n ) )./stdcell.LAV;
                    f_ext( vert(1) ) = f_ext( vert(1) ) + fval;
                    fval = area./2 .* ( stdcell.wq' ...
                        * FunHandle( rq2, n ) )./stdcell.LAV;
                    f_ext( vert(2) ) = f_ext( vert(2) ) + fval;
                end
                f = ( stdcell.P * ( FunHandle( r, n ) ) ) .*stdcell.LAVToCV;
                test.verifyEqual( f, f_ext, 'AbsTol', test.tol);
            end
        end% func
        
        function TestTriProjectMatrix( test )
            stdcell = test.cell;
            r = stdcell.r;
            s = stdcell.s;
            FunHandle = @( r,s,n ) ( r - s ).^(n);
            for n = 0:stdcell.N
                f_ext = zeros( stdcell.Np, 1 );
                for i = 1:stdcell.NFV
                    vert = stdcell.EToV(:, i);
                    rv = r( vert );
                    sv = s( vert );
                    rc = mean( rv ); 
                    sc = mean( sv );
                    area = TriArea( rv, sv );
                    for f = 1:stdcell.Nface
                        n1 = f;
                        n2 = mod(f, stdcell.Nface)+1;
                        v1 = vert(n1);
                        v2 = vert(n2);
                        r1 = r(v1); r2 = r(v2);
                        s1 = s(v1); s2 = s(v2);
                        rf = ( r1 + r2 )*0.5;
                        sf = ( s1 + s2 )*0.5;
                        rv1 = [rc, r1, rf]';
                        sv1 = [sc, s1, sf]';
                        rq1 = stdcell.project_vert2quad(rv1);
                        sq1 = stdcell.project_vert2quad(sv1);
                        fval = area./6 .* ( stdcell.wq' ...
                            * FunHandle( rq1, sq1, n ) )./stdcell.LAV;
                        f_ext( v1 ) = f_ext( v1 ) + fval;
                        
                        rv2 = [rc, rf, r2]';
                        sv2 = [sc, sf, s2]';
                        rq2 = stdcell.project_vert2quad(rv2);
                        sq2 = stdcell.project_vert2quad(sv2);
                        fval = area./6.*( stdcell.wq' ...
                            * FunHandle( rq2, sq2, n ) )./stdcell.LAV;
                        f_ext( v2 ) = f_ext( v2 ) + fval;
                    end
                end
                
                f = ( stdcell.P * ( FunHandle( r, s, n ) ) ) .*stdcell.LAVToCV;
                test.verifyEqual( f, f_ext, 'AbsTol', test.tol);
            end
        end% func
    end
    
end

function area = TriArea(vx, vy)
% calculate the area of triangle from given vertex coordinate
a = sqrt( (vx(1,:) - vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2 );
b = sqrt( (vx(2,:) - vx(3,:)).^2 + (vy(2,:)-vy(3,:)).^2 );
c = sqrt( (vx(3,:) - vx(1,:)).^2 + (vy(3,:)-vy(1,:)).^2 );
p = (a+b+c)./2;
area = sqrt(p.*(p-a).*(p-b).*(p-c));
end

