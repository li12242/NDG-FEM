%> @brief Unit test for the StdCell classes.
%
%> @code
%>     addpath(pwd);
%>     results = runtests('NdgCell/test/StdCellTest.m');
%> @endcode
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef StdCellTest < matlab.unittest.TestCase
    
    properties(MethodSetupParameter)
        %> test cell types
        type = {...
%             StdCellType.Line, ...
            StdCellType.Tri, ...
%             StdCellType.Quad
            }
        %> test cell orders
        order = {3,4,5}
    end
    
    properties(Constant)
        %> tolerance
        tol = 1e-9;
    end
    
    properties
        % cell object
        cell
    end
    
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_std_cell(test, type, order)
            test.cell = getStdCell(order, type);
        end% func
    end
    
    methods(Test, ParameterCombination = 'sequential')
        function test_point_coor(test)
            [r_ext, s_ext, t_ext] = get_ext_coor(test.cell);
            test.verifyEqual(test.cell.r, r_ext, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.s, s_ext, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.t, t_ext, 'AbsTol', test.tol);
        end
        
        function test_vand_matrix(test)
            [ vand_ext ] = get_ext_vandmatrix(test.cell);
            test.verifyEqual(test.cell.V, vand_ext, 'AbsTol', test.tol);
        end
        
        function test_deri_matrix(test)
            [ dr, ds, dt ] = get_ext_derimatrix(test.cell);
            test.verifyEqual(test.cell.Dr, dr, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.Ds, ds, 'AbsTol', test.tol);
            test.verifyEqual(test.cell.Dt, dt, 'AbsTol', test.tol);
        end
        
        function test_quadrature_weight(test)
            N = test.cell.N;
            r = test.cell.rq;
            s = test.cell.sq;
            t = test.cell.tq;
            w = test.cell.wq;
            switch test.cell.type
                case ndg_lib.std_cell_type.Line
                    for i = 1:(N*2-1)
                        f = r.^i;
                        int_val = sum( w.*f );
                        % the exact value for volume integral
                        ext_val = (1 - (-1)^(i+1))/( i+1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case ndg_lib.std_cell_type.Tri
                    for i = 1:(N*2-1)
                        f = r.^i;
                        int_val = sum( w.*f );
                        % the exact value for volume integral
                        ext_val = (-1).^i * 2/( i + mod(i,2) + 1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case ndg_lib.std_cell_type.Quad
            end
        end% func
        
        %> test the orthgonality of the basis function
        function test_orthgonal_func(test)
            N = test.cell.N;
            w = test.cell.wq;
            V = zeros(test.cell.Nq, test.cell.Np);
            for n = 1:test.cell.Np
                V(:, n) = test.cell.orthogonal_func(...
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
        
        function test_derivative_matrix(test)
            Dr = test.cell.Dr;
            Ds = test.cell.Ds;
            Dt = test.cell.Dt;
            
            r = test.cell.r;
            s = test.cell.s;
            t = test.cell.t;
            switch test.cell.type
                case StdCellType.Line
                    test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                case {StdCellType.Tri, StdCellType.Quad}
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
        
    end% methods
end% classdef

function [Dr, Ds, Dt] = get_ext_derimatrix(cell)
Dr = zeros(cell.Np, cell.Np);
Ds = zeros(cell.Np, cell.Np);
Dt = zeros(cell.Np, cell.Np);
switch cell.type
    case StdCellType.Line
        folder = 'line/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
    case StdCellType.Tri
        folder = 'tri/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
        Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
    case StdCellType.Quad
        folder = 'quad/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
        Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
end
end% func

function [r_ext, s_ext, t_ext] = get_ext_coor(cell)
r_ext = zeros(cell.Np, 1);
s_ext = zeros(cell.Np, 1);
t_ext = zeros(cell.Np, 1);

switch cell.type
    case StdCellType.Line
        folder = 'line/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
    case StdCellType.Tri
        folder = 'tri/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
    case StdCellType.Quad
        folder = 'quad/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
end
end% func

function [ V_ext ] = get_ext_vandmatrix(cell)
switch cell.type
    case StdCellType.Line
        folder = 'line/Vand_Test/';
    case StdCellType.Tri
        folder = 'tri/Vand_Test/';
    case StdCellType.Quad
        folder = 'quad/Vand_Test/';
end
V_ext = load([folder, 'Vand_', num2str(cell.N), '.cc']);
end
