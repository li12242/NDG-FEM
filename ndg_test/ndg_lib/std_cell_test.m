classdef std_cell_test < matlab.unittest.TestCase  
    %STD_CELL_TEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(MethodSetupParameter)
        type = {ndg_lib.std_cell_type.Line, ...
            ndg_lib.std_cell_type.Tri, ...
            ndg_lib.std_cell_type.Quad
            }
        order = {3,4,5}
    end
    
    properties(Constant)
        tol = 1e-9;
    end
    
    properties
        cell
    end
    
    methods(TestMethodSetup)
        function set_std_cell(test, type, order) % setup the std_cell obj
            test.cell = ndg_lib.get_std_cell(order, type);
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
            r = test.cell.r;
            s = test.cell.s;
            t = test.cell.t;
            w = test.cell.w;
            switch test.cell.type
                case ndg_lib.std_cell_type.Line
                case ndg_lib.std_cell_type.Tri
                    for i = 1:test.cell.N
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
        
        function test_orthgonal_func(test)
        % test the orthgonality of the basis function
            w = test.cell.w;
            V = test.cell.V;
            
            switch test.cell.type
                case ndg_lib.std_cell_type.Tri
                    % the maximum exact integral degree of polynomial
                    maxDeg = test.cell.N + 2; 
                otherwise
                    maxDeg = test.cell.N*2 - 1;
            end% switch
            
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
                case ndg_lib.std_cell_type.Line
                    test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
                        'AbsTol', test.tol );
                case {ndg_lib.std_cell_type.Tri, ...
                        ndg_lib.std_cell_type.Quad}
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
    case ndg_lib.std_cell_type.Line
        folder = 'std_cell/line/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
    case ndg_lib.std_cell_type.Tri
        folder = 'std_cell/tri/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
        Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
    case ndg_lib.std_cell_type.Quad  
        folder = 'std_cell/quad/DrDs_Test/';
        Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
        Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
end
end% func

function [r_ext, s_ext, t_ext] = get_ext_coor(cell)
r_ext = zeros(cell.Np, 1);
s_ext = zeros(cell.Np, 1);
t_ext = zeros(cell.Np, 1);

switch cell.type
    case ndg_lib.std_cell_type.Line
        folder = 'std_cell/line/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
    case ndg_lib.std_cell_type.Tri
        folder = 'std_cell/tri/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
    case ndg_lib.std_cell_type.Quad  
        folder = 'std_cell/quad/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
end
end% func

function [ V_ext ] = get_ext_vandmatrix(cell)
switch cell.type
    case ndg_lib.std_cell_type.Line
        folder = 'std_cell/line/Vand_Test/';
    case ndg_lib.std_cell_type.Tri
        folder = 'std_cell/tri/Vand_Test/';
    case ndg_lib.std_cell_type.Quad  
        folder = 'std_cell/quad/Vand_Test/';
end
V_ext = load([folder, 'Vand_', num2str(cell.N), '.cc']);
end
