%> @brief Unit tests for the pure-Matlab Polylib gateway ports.
%>
%> Two groups of checks are performed:
%>  - property tests: orthonormalization, known Legendre values, finite
%>    difference derivatives, LGL symmetry/exactness (independent of mex);
%>  - golden regression: comparison against the mex baseline captured by
%>    genGoldenPolylib.m BEFORE the mex binaries were deleted
%>    (golden cases are limited to n<=12, see genGoldenPolylib).
%>
% ======================================================================
%> This class is part of the NDG-FEM software.
% ======================================================================
classdef PolylibGatewayTest < matlab.unittest.TestCase

    properties(Constant)
        %> tolerance
        tol = 1e-12;
    end

    properties
        %> mex baseline data
        golden
    end

    methods(TestClassSetup)
        function addPathAndLoadGolden(testCase)
            here = fileparts(mfilename('fullpath'));
            addpath( fullfile(here, '..', '..', 'lib') );
            data = load( fullfile(here, 'PolylibGoldenTest', ...
                'polylib_golden.mat'), 'golden' );
            testCase.golden = data.golden;
        end
    end

    methods(Test)
        %> orthonormality: \int_{-1}^{1} \tilde P_n^2 dr = 1
        function testJacobiPNormalization(testCase)
            [x, w] = zwglj(30);
            for n = 0:8
                P = JacobiP(x, 0, 0, n);
                testCase.verifyEqual( sum(w.*P.*P), 1, ...
                    'AbsTol', testCase.tol );
            end
        end

        %> known values of the normalized Legendre polynomial
        function testJacobiPKnownLegendre(testCase)
            x = linspace(-1, 1, 11)';
            testCase.verifyEqual( JacobiP(x, 0, 0, 0), ...
                1/sqrt(2)*ones(11, 1), 'AbsTol', testCase.tol );
            testCase.verifyEqual( JacobiP(x, 0, 0, 1), ...
                sqrt(1.5)*x, 'AbsTol', testCase.tol );
            testCase.verifyEqual( JacobiP(x, 0, 0, 2), ...
                sqrt(2.5)*0.5*(3*x.^2 - 1), 'AbsTol', testCase.tol );
        end

        %> the mex always returns a numel(r)-by-1 column
        function testJacobiPColumnOutput(testCase)
            r = linspace(-1, 1, 7); % row vector on purpose
            P = JacobiP(r, 0, 0, 3);
            testCase.verifySize(P, [7, 1]);
        end

        function testGradJacobiPZeroOrder(testCase)
            testCase.verifyEqual( GradJacobiP(linspace(-1, 1, 5)', 0, 0, 0), ...
                zeros(5, 1) );
        end

        %> finite-difference check of the derivative identity
        function testGradJacobiPFiniteDiff(testCase)
            h = 1e-6;
            x = linspace(-0.9, 0.9, 9)';
            for ab = [0, 0; 2, 0; 3, 1]'
                for n = 1:5
                    fd = ( JacobiP(x + h, ab(1), ab(2), n) ...
                         - JacobiP(x - h, ab(1), ab(2), n) )/(2*h);
                    testCase.verifyEqual( GradJacobiP(x, ab(1), ab(2), n), ...
                        fd, 'AbsTol', 1e-5, 'RelTol', 1e-5 );
                end
            end
        end

        function testZwgljDegenerateCase(testCase)
            [z, w] = zwglj(1);
            testCase.verifyEqual(z, 0);
            testCase.verifyEqual(w, 2);
        end

        %> LGL properties: endpoints, ascending, symmetry, weight sums
        function testZwgljProperty(testCase)
            for np = 2:8
                [z, w] = zwglj(np);
                testCase.verifyEqual(z(1), -1);
                testCase.verifyEqual(z(end), 1);
                testCase.verifyTrue( issorted(z) );
                testCase.verifyEqual(z, -flipud(z), 'AbsTol', testCase.tol);
                testCase.verifyEqual(w, flipud(w), 'AbsTol', testCase.tol);
                testCase.verifyEqual(sum(w), 2, 'AbsTol', testCase.tol);
            end
        end

        %> the np-point LGL rule integrates polynomials of degree 2*np-3 exactly
        function testZwgljExactness(testCase)
            np = 8;
            [z, w] = zwglj(np);
            for k = 0:(2*np - 3)
                exact = (1 - (-1)^(k + 1))/(k + 1);
                testCase.verifyEqual( w.'*(z.^k), exact, ...
                    'AbsTol', testCase.tol );
            end
        end

        % ---------------- golden regression (mex baseline) ----------------

        function testJacobiPGolden(testCase)
            for c = 1:numel(testCase.golden.jacobiP)
                ref = testCase.golden.jacobiP(c);
                testCase.verifyEqual( ...
                    JacobiP(ref.r, ref.alpha, ref.beta, ref.n), ref.P, ...
                    'AbsTol', testCase.tol, 'RelTol', testCase.tol, ...
                    sprintf('JacobiP golden mismatch: alpha=%g beta=%g n=%d', ...
                    ref.alpha, ref.beta, ref.n) );
            end
        end

        function testGradJacobiPGolden(testCase)
            for c = 1:numel(testCase.golden.gradJacobiP)
                ref = testCase.golden.gradJacobiP(c);
                testCase.verifyEqual( ...
                    GradJacobiP(ref.r, ref.alpha, ref.beta, ref.n), ref.dP, ...
                    'AbsTol', testCase.tol, 'RelTol', testCase.tol, ...
                    sprintf('GradJacobiP golden mismatch: alpha=%g beta=%g n=%d', ...
                    ref.alpha, ref.beta, ref.n) );
            end
        end

        function testZwgljGolden(testCase)
            for c = 1:numel(testCase.golden.zwglj)
                ref = testCase.golden.zwglj(c);
                [z, w] = zwglj(ref.np);
                testCase.verifyEqual(z, ref.z, 'AbsTol', testCase.tol, ...
                    sprintf('zwglj z mismatch: np=%d', ref.np) );
                testCase.verifyEqual(w, ref.w, 'AbsTol', testCase.tol, ...
                    'RelTol', testCase.tol, ...
                    sprintf('zwglj w mismatch: np=%d', ref.np) );
            end
        end
    end
end
