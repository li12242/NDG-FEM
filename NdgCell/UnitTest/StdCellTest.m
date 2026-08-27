%> @brief Unit test for the StdCell classes.
%>
%> Covers StdLine/StdTri/StdQuad at orders 3-5:
%>  - interpolation node coordinates vs reference data (*.cc);
%>  - Vandermonde matrix vs reference data (*.cc);
%>  - derivative matrix Dr/Ds vs reference data (*.cc);
%>  - derivative identity (e.g. Dr*r = 1);
%>  - quadrature weights against exact integrals;
%>  - orthonormality of the basis w.r.t. the quadrature;
%>  - face node indexing (Fmask);
%>  - consistency of face-cell quadrature interpolation.
%>
% ======================================================================
%> This class is part of the NDG-FEM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef StdCellTest < matlab.unittest.TestCase

    properties(MethodSetupParameter)
        %> test cell types
        type = {...
            enumStdCell.Line, ...
            enumStdCell.Tri, ...
            enumStdCell.Quad
            }
        %> test cell orders
        order = { 3,4,5 }
    end

    properties(Constant)
        %> tolerance
        tol = 1e-9;
    end

    properties
        %> folder of this test file (reference data lives below it)
        testdir
        %> cell object
        cell
    end

    methods(TestClassSetup)
        %> make the NdgCell module and the Polylib ports reachable
        function addModulePath(testCase)
            here = fileparts(mfilename('fullpath'));
            testCase.testdir = here;
            addpath( fullfile(here, '..') );              % NdgCell
            addpath( fullfile(here, '..', '..', 'lib') ); % lib
        end
    end

    methods(TestMethodSetup)
        %> get the StdCell object
        function setStdCell(test, type, order)
            test.cell = getStdCell(order, type);
        end% func
    end

    methods(Test, ParameterCombination = 'sequential')
        function testPointCoor(test)
            [ r_ext, s_ext, t_ext ] = getExtCoor(test.testdir, test.cell);
            test.verifyEqual(test.cell.r, r_ext, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.s, s_ext, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.t, t_ext, 'AbsTol', test.tol);
        end

        %> @brief test the basis function values on edge quadrature points.
        %> Test whether the values of the basis function on IPPS at edges
        %> are equal to the basis function values from the facial cell Vq.
        %> \f$ l_i^{\Omega}(r_q, s_q) = l_{fi}^{\partial \Omega}(r_q), \f$
        %> while \f$ l_i^{\Omega} \f$ is the nodal basis function of
        %> interpolation nodes \f$(r_i, s_i)\f$, and
        %> \f$ \l_{fi}^{\partial \Omega} \f$ are the nodal basis function
        %> of the boundary cell, and on the interpolation nodel \f$(r_{fi})\f$.
        function testEdgeBasisVal(test)
            facetype = test.cell.faceType;
            for f = 1:test.cell.Nface
                bcell = getStdCell(test.cell.N, facetype(f));

                r = test.cell.r(test.cell.Fmask(:,f));
                s = test.cell.s(test.cell.Fmask(:,f));
                rq = bcell.projectNode2Quad(r);
                sq = bcell.projectNode2Quad(s);
                Vq = zeros(bcell.Nq, test.cell.Np);
                ind = test.cell.Fmask(:, f);
                for n = 1:test.cell.Np
                    Vq(ind, n) = test.cell.evaluateOrthogonalFunc(bcell.N, n, rq, sq);
                end% func
                Vq = Vq / (test.cell.V);
                test.verifyEqual(Vq(ind, ind), bcell.Vq, 'AbsTol', test.tol);
            end
        end

        function testVandMatrix(test)
            [ vand_ext ] = getExtVandMatrix(test.testdir, test.cell);
            test.verifyEqual(test.cell.V, vand_ext, 'AbsTol', test.tol);
        end

        %> derivative matrices vs reference data (*.cc)
        function testDeriMatrix(test)
            [ dr, ds, dt ] = getExtDeriMatrix(test.testdir, test.cell);
            test.verifyEqual(test.cell.Dr, dr, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.Ds, ds, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.Dt, dt, 'AbsTol', test.tol);
        end

        %> derivative identity: d/dr of the linear function r is exactly 1
        function testDerivativeIdentity(test)
            Dr = test.cell.Dr;
            Ds = test.cell.Ds;
            r = test.cell.r;
            s = test.cell.s;
            switch test.cell.type
                case enumStdCell.Line
                    test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                case {enumStdCell.Tri, enumStdCell.Quad}
                    test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                    test.verifyEqual( Ds*s, ones(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                    test.verifyEqual( Dr*s, zeros(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                    test.verifyEqual( Ds*r, zeros(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
            end
        end% func

        %> quadrature weights against exact monomial integrals
        function testQuadratureWeight(test)
            N = test.cell.N;
            r = test.cell.rq;
            w = test.cell.wq;
            switch test.cell.type
                case enumStdCell.Line
                    for i = 1:(N*2-1)
                        int_val = sum( w.*(r.^i) );
                        % the exact value for volume integral
                        ext_val = (1 - (-1)^(i+1))/( i+1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case enumStdCell.Tri
                    for i = 1:(N*2-1)
                        int_val = sum( w.*(r.^i) );
                        % the exact value for volume integral
                        ext_val = (-1)^i * 2/( i + mod(i,2) + 1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case enumStdCell.Quad
                    %> \iint_{[-1,1]^2} r^i dr ds = 2(1-(-1)^{i+1})/(i+1)
                    for i = 1:(N*2-1)
                        int_val = sum( w.*(r.^i) );
                        ext_val = 2*(1 - (-1)^(i+1))/( i+1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
            end
        end% func

        %> test the orthgonality of the basis function
        function testOrthogonalFunc(test)
            N = test.cell.N;
            w = test.cell.wq;
            V = zeros(test.cell.Nq, test.cell.Np);
            for n = 1:test.cell.Np
                V(:, n) = test.cell.evaluateOrthogonalFunc(...
                    N, n, test.cell.rq, test.cell.sq, test.cell.tq);
            end

            % the maximum exact integral degree of polynomial
            maxDeg = test.cell.N*2 - 1;

            for i = 1:test.cell.Np
                for j = 1:test.cell.Np
                    if( (i+j)>maxDeg )
                        continue;
                    end
                    temp = sum( w.*V(:,i).*V(:,j) );
                    test.verifyEqual( temp, double(i==j), ...
                        'AbsTol', test.tol );
                end
            end
        end% func

        %> face node index: valid range and counts per face
        function testFmask(test)
            for f = 1:test.cell.Nface
                ind = test.cell.Fmask(:, f);
                nNode = nnz(ind);
                test.verifyEqual( nNode, test.cell.Nfp(f) );
                test.verifyTrue( all( ind(ind>0) >= 1 & ...
                    ind(ind>0) <= test.cell.Np ) );
            end
        end

    end% methods
end% classdef

function [r_ext, s_ext, t_ext] = getExtCoor(testdir, cell)
r_ext = zeros(cell.Np, 1);
s_ext = zeros(cell.Np, 1);
t_ext = zeros(cell.Np, 1);

switch cell.type
    case enumStdCell.Line
        folder = fullfile(testdir, 'StdLineTest', 'Coor_Test');
        r_ext = load( fullfile(folder, ['r_', num2str(cell.N), '.cc']) );
    case enumStdCell.Tri
        folder = fullfile(testdir, 'StdTriTest', 'Coor_Test');
        r_ext = load( fullfile(folder, ['r_', num2str(cell.N), '.cc']) );
        s_ext = load( fullfile(folder, ['s_', num2str(cell.N), '.cc']) );
    case enumStdCell.Quad
        folder = fullfile(testdir, 'StdQuadTest', 'Coor_Test');
        r_ext = load( fullfile(folder, ['r_', num2str(cell.N), '.cc']) );
        s_ext = load( fullfile(folder, ['s_', num2str(cell.N), '.cc']) );
end
end% func

function [ V_ext ] = getExtVandMatrix(testdir, cell)
switch cell.type
    case enumStdCell.Line
        folder = fullfile(testdir, 'StdLineTest', 'Vand_Test');
    case enumStdCell.Tri
        folder = fullfile(testdir, 'StdTriTest', 'Vand_Test');
    case enumStdCell.Quad
        folder = fullfile(testdir, 'StdQuadTest', 'Vand_Test');
end
V_ext = load( fullfile(folder, ['Vand_', num2str(cell.N), '.cc']) );
end% func

function [Dr, Ds, Dt] = getExtDeriMatrix(testdir, cell)
Dr = zeros(cell.Np, cell.Np);
Ds = zeros(cell.Np, cell.Np);
Dt = zeros(cell.Np, cell.Np);
switch cell.type
    case enumStdCell.Line
        folder = fullfile(testdir, 'StdLineTest', 'DrDs_Test');
        Dr = load( fullfile(folder, ['Dr_', num2str(cell.N), '.cc']) );
    case enumStdCell.Tri
        folder = fullfile(testdir, 'StdTriTest', 'DrDs_Test');
        Dr = load( fullfile(folder, ['Dr_', num2str(cell.N), '.cc']) );
        Ds = load( fullfile(folder, ['Ds_', num2str(cell.N), '.cc']) );
    case enumStdCell.Quad
        folder = fullfile(testdir, 'StdQuadTest', 'DrDs_Test');
        Dr = load( fullfile(folder, ['Dr_', num2str(cell.N), '.cc']) );
        Ds = load( fullfile(folder, ['Ds_', num2str(cell.N), '.cc']) );
end
end% func
