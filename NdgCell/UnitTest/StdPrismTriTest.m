%> @brief Unit test for the StdPrismTri class.
%>
%> StdPrismTri is constructed directly (StdPrismTri(Nh, Nz)) with equal
%> horizontal/vertical orders here; the checks are property based (no
%> external reference data):
%>  - node counts and layout (tensor of triangle and 1D LGL nodes);
%>  - face node counts and index validity;
%>  - quadrature weights against exact monomial integrals;
%>  - orthonormality of low-order modes w.r.t. the quadrature;
%>  - derivative identities (Dr*r = 1 etc.);
%>  - vertex projection of a linear field.
%>
% ======================================================================
%> This class is part of the NDG-FEM software.
% ======================================================================
classdef StdPrismTriTest < matlab.unittest.TestCase

    properties(MethodSetupParameter)
        %> test cell orders (Nh = Nz)
        order = { 2,3,4 }
    end

    properties(Constant)
        %> tolerance
        tol = 1e-9;
    end

    properties
        %> cell object
        cell
    end

    methods(TestClassSetup)
        function addModulePath(testCase)
            here = fileparts(mfilename('fullpath'));
            addpath( fullfile(here, '..') );              % NdgCell
            addpath( fullfile(here, '..', '..', 'lib') ); % lib
        end
    end

    methods(TestMethodSetup)
        function set_prism(test, order)
            test.cell = StdPrismTri(order, order);
        end
    end

    methods(Test, ParameterCombination = 'sequential')
        %> node counts and tensor layout
        function test_node_layout(test)
            N = test.cell.N;
            test.verifyEqual( test.cell.Nph, (N+1)*(N+2)/2 );
            test.verifyEqual( test.cell.Npz, N+1 );
            test.verifyEqual( test.cell.Np, test.cell.Nph*test.cell.Npz );
            test.verifySize( test.cell.r, [test.cell.Np, 1] );
            test.verifySize( test.cell.s, [test.cell.Np, 1] );
            test.verifySize( test.cell.t, [test.cell.Np, 1] );
            %> horizontal nodes repeat per vertical level
            test.verifyEqual( test.cell.r(1:test.cell.Nph), ...
                test.cell.r(test.cell.Nph+1:2*test.cell.Nph), ...
                'AbsTol', test.tol );
            %> vertical levels are sorted LGL values: first level at t = -1
            test.verifyEqual( test.cell.t(1:test.cell.Nph), ...
                -ones(test.cell.Nph, 1), 'AbsTol', test.tol );
        end

        %> face node counts and index validity
        function test_fmask(test)
            N = test.cell.N;
            for f = 1:3 % quad side faces
                test.verifyEqual( test.cell.Nfp(f), (N+1)*(N+1) );
            end
            for f = 4:5 % tri top/bottom faces
                test.verifyEqual( test.cell.Nfp(f), test.cell.Nph );
            end
            test.verifyEqual( test.cell.TNfp, sum(test.cell.Nfp) );
            for f = 1:test.cell.Nface
                ind = test.cell.Fmask(:, f);
                test.verifyEqual( nnz(ind), test.cell.Nfp(f) );
                test.verifyTrue( all( ind(ind>0) >= 1 & ...
                    ind(ind>0) <= test.cell.Np ) );
            end
        end

        %> quadrature weights against exact monomial integrals
        %> \iiint r dV = -4/3 (centroid x = -1/3, Area(tri) = 2, height = 2),
        %> \iiint t dV = 0, total volume = 4.
        function test_quadrature_weight(test)
            w = test.cell.wq;
            test.verifyEqual( sum(w), 4, 'AbsTol', test.tol );
            test.verifyEqual( sum( w.*test.cell.rq ), -4/3, ...
                'AbsTol', test.tol );
            test.verifyEqual( sum( w.*test.cell.sq ), -4/3, ...
                'AbsTol', test.tol );
            test.verifyEqual( sum( w.*test.cell.tq ), 0, ...
                'AbsTol', test.tol );
        end

        %> orthonormality of modes with combined degree <= N-1
        %> (triquad and LGL rules are exact up to degree 2N-1)
        function test_orthgonal_func(test)
            N = test.cell.N;
            w = test.cell.wq;
            rq = test.cell.rq; sq = test.cell.sq; tq = test.cell.tq;
            [ a, b ] = ndgcell.rstoab( rq, sq );
            Nph = test.cell.Nph;
            modal = zeros(test.cell.Nq, test.cell.Np);
            deg = zeros(test.cell.Np, 1);
            for n = 1:test.cell.Np
                td1 = mod( n-1, Nph ) + 1;
                td2 = ceil( n/Nph );
                [ i1, j1 ] = ndgcell.trans_ind( N, td1 );
                modal(:, n) = ndgcell.simplex2DP( a, b, i1, j1 ) ...
                    .* JacobiP( tq, 0, 0, td2-1 );
                deg(n) = i1 + j1 + (td2 - 1);
            end
            sel = find( deg <= N-1 );
            for ii = 1:numel(sel)
                for jj = ii:numel(sel)
                    temp = sum( w.*modal(:,sel(ii)).*modal(:,sel(jj)) );
                    test.verifyEqual( temp, double(sel(ii)==sel(jj)), ...
                        'AbsTol', test.tol );
                end
            end
        end

        %> derivative identities of linear functions
        function test_derivative_identity(test)
            Np = test.cell.Np;
            test.verifyEqual( test.cell.Dr*test.cell.r, ones(Np, 1), ...
                'AbsTol', test.tol );
            test.verifyEqual( test.cell.Ds*test.cell.s, ones(Np, 1), ...
                'AbsTol', test.tol );
            test.verifyEqual( test.cell.Dt*test.cell.t, ones(Np, 1), ...
                'AbsTol', test.tol );
            test.verifyEqual( test.cell.Dr*test.cell.s, zeros(Np, 1), ...
                'AbsTol', test.tol );
            test.verifyEqual( test.cell.Ds*test.cell.t, zeros(Np, 1), ...
                'AbsTol', test.tol );
        end

        %> vertex projection is exact for a linear field
        function test_project_vert2node(test)
            vert_val = test.cell.vr + test.cell.vs + test.cell.vt; % 6-by-1
            node_val = test.cell.project_vert2node(vert_val);
            ext_val = test.cell.r + test.cell.s + test.cell.t;
            test.verifyEqual( node_val, ext_val, 'AbsTol', test.tol );
        end
    end
end
