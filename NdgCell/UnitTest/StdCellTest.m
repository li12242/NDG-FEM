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
            NdgCellType.Line, ...
            NdgCellType.Tri, ...
            NdgCellType.Quad
            }
        %> test cell orders
        order = {3,4,5}
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
        
        %> @brief test the basis function values on edge quadrature points.
        %> Test whether the values of the basis function on IPPS at edges  
        %> are equal to the basis function values from the facial cell Vq.
        %> \f$ l_i^{\Omega}(r_q, s_q) = l_{fi}^{\partial \Omega}(r_q), \f$
        %> while \f$ l_i^{\Omega} \f$ is the nodal basis function of
        %> interpolation nodes \f$(r_i, s_i)\f$, and 
        %> \f$ \l_{fi}^{\partial \Omega} \f$ are the nodal basis function
        %> of the boundary cell, and on the interpolation nodel \f$(r_{fi})\f$.
        function test_edge_basis_val(test)
            facetype = test.cell.faceType;
            for f = 1:test.cell.Nface
                bcell = getStdCell(test.cell.N, facetype(f));
                
                r = test.cell.r(test.cell.Fmask(:,f));
                s = test.cell.s(test.cell.Fmask(:,f));
                rq = bcell.project_node2quad(r);
                sq = bcell.project_node2quad(s);
                Vq = zeros(bcell.Nq, test.cell.Np);
                ind = test.cell.Fmask(:, f);
                for n = 1:test.cell.Np
                    Vq(ind, n) = test.cell.orthogonal_func(bcell.N, n, rq, sq);
                end% func
                Vq = Vq / (test.cell.V);
                bcell = getStdCell(test.cell.N, facetype(f));
                test.verifyEqual(Vq(ind, ind), bcell.Vq, 'AbsTol', test.tol);
            end
        end
        
        function test_vand_matrix(test)
            [ vand_ext ] = get_ext_vandmatrix(test.cell);
            test.verifyEqual(test.cell.V, vand_ext, 'AbsTol', test.tol);
        end
        
%         function test_deri_matrix(test)
%             [ dr, ds, dt ] = get_ext_derimatrix(test.cell);
%             test.verifyEqual(test.cell.Dr, dr, 'AbsTol', test.tol);
%             test.verifyEqual(test.cell.Ds, ds, 'AbsTol', test.tol);
%             test.verifyEqual(test.cell.Dt, dt, 'AbsTol', test.tol);
%         end
        
        function test_quadrature_weight(test)
            N = test.cell.N;
            r = test.cell.rq;
            s = test.cell.sq;
            t = test.cell.tq;
            w = test.cell.wq;
            switch test.cell.type
                case NdgCellType.Line
                    for i = 1:(N*2-1)
                        f = r.^i;
                        int_val = sum( w.*f );
                        % the exact value for volume integral
                        ext_val = (1 - (-1)^(i+1))/( i+1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case NdgCellType.Tri
                    for i = 1:(N*2-1)
                        f = r.^i;
                        int_val = sum( w.*f );
                        % the exact value for volume integral
                        ext_val = (-1).^i * 2/( i + mod(i,2) + 1 );
                        test.verifyEqual( ...
                            int_val, ext_val, 'AbsTol', test.tol );
                    end
                case NdgCellType.Quad
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
        
%         function test_derivative_matrix(test)
%             Dr = test.cell.Dr;
%             Ds = test.cell.Ds;
%             Dt = test.cell.Dt;
%             
%             r = test.cell.r;
%             s = test.cell.s;
%             t = test.cell.t;
%             switch test.cell.type
%                 case NdgCellType.Line
%                     test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
%                         'AbsTol', test.tol );
%                 case {NdgCellType.Tri, NdgCellType.Quad}
%                     test.verifyEqual( Dr*r, ones(test.cell.Np, 1), ...
%                         'AbsTol', test.tol );
%                     test.verifyEqual( Ds*s, ones(test.cell.Np, 1), ...
%                         'AbsTol', test.tol );
%                     test.verifyEqual( Dr*s, zeros(test.cell.Np, 1), ...
%                         'AbsTol', test.tol );
%                     test.verifyEqual( Ds*r, zeros(test.cell.Np, 1), ...
%                         'AbsTol', test.tol );
%             end
%         end% func
        
    end% methods
end% classdef

% function [Dr, Ds, Dt] = get_ext_derimatrix(cell)
% Dr = zeros(cell.Np, cell.Np);
% Ds = zeros(cell.Np, cell.Np);
% Dt = zeros(cell.Np, cell.Np);
% switch cell.type
%     case NdgCellType.Line
%         folder = 'StdLineTest/DrDs_Test/';
%         Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
%     case NdgCellType.Tri
%         folder = 'StdTriTest/DrDs_Test/';
%         Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
%         Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
%     case NdgCellType.Quad
%         folder = 'StdQuadTest/DrDs_Test/';
%         Dr = load([folder, 'Dr_', num2str(cell.N), '.cc']);
%         Ds = load([folder, 'Ds_', num2str(cell.N), '.cc']);
% end
% end% func

function [r_ext, s_ext, t_ext] = get_ext_coor(cell)
r_ext = zeros(cell.Np, 1);
s_ext = zeros(cell.Np, 1);
t_ext = zeros(cell.Np, 1);

switch cell.type
    case NdgCellType.Line
        folder = 'StdLineTest/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
    case NdgCellType.Tri
        folder = 'StdTriTest/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
    case NdgCellType.Quad
        folder = 'StdQuadTest/Coor_Test/';
        r_ext = load([folder, 'r_', num2str(cell.N), '.cc']);
        s_ext = load([folder, 's_', num2str(cell.N), '.cc']);
end
end% func

function [ V_ext ] = get_ext_vandmatrix(cell)
switch cell.type
    case NdgCellType.Line
        folder = 'StdLineTest/Vand_Test/';
    case NdgCellType.Tri
        folder = 'StdTriTest/Vand_Test/';
    case NdgCellType.Quad
        folder = 'StdQuadTest/Vand_Test/';
end
V_ext = load([folder, 'Vand_', num2str(cell.N), '.cc']);
end
