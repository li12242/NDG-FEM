%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgMeshTest < matlab.unittest.TestCase
    %UMESHTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(MethodSetupParameter)
        type = {...
            %NdgCellType.Line, ...
            NdgCellType.Tri, ...
            NdgCellType.Quad, ...
            }
        order = {1, 2, 3}
    end
    
    properties
        cell
        mesh
        tol = 1e-9;
    end
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_std_cell(test, type, order)
            test.cell = getStdCell(order, type);
        end% func
        
        function set_test_cell(test)
            stdCell = test.cell;
            switch test.cell.type
                case NdgCellType.Line
                    Nv = stdCell.Nv;
                    vx = stdCell.vr./2+1;
                    K = 1;
                    EToV = [1,2]';
                    EToR = 1;
                    BCToV = [1, 0; 2, 0]';
                    testMesh = NdgMesh1d(stdCell,Nv,vx,K,EToV,EToR,BCToV);
                case NdgCellType.Tri
                    Nv = stdCell.Nv;
                    vx = stdCell.vr/2+3;
                    vy = stdCell.vs/2+4;
                    K = 1;
                    EToV = [1,2,3]';
                    EToR = 1;
                    BCToV = [1,2,0; 2,3,0; 3,1,0]';
                    testMesh = NdgMesh2d(stdCell,Nv,vx,vy,K,EToV,EToR,BCToV);
                case NdgCellType.Quad 
                    Nv = stdCell.Nv;
                    vx = stdCell.vr*3 + 2*stdCell.vs + 2;
                    vy = stdCell.vs/2+1;
                    K = 1;
                    EToV = [1,2,3,4]';
                    EToR = 1;
                    BCToV = [1,2,0; 2,3,0; 3,4,0; 4,1,0]';
                    testMesh = NdgMesh2d(stdCell,Nv,vx,vy,K,EToV,EToR,BCToV);
            end
            test.mesh = testMesh;
        end% func
    end% methods
    
    methods(Test, ParameterCombination = 'sequential')
        function test_coor(test)
            r = test.cell.r;
            s = test.cell.s;
            t = test.cell.t;
            
            x = test.mesh.x;
            y = test.mesh.y;
            z = test.mesh.z;
            switch test.cell.type
                case NdgCellType.Line
                    test.verifyEqual( r/2 + 1, x, 'AbsTol', test.tol);
                case NdgCellType.Tri
                    test.verifyEqual( r/2 + 3, x, 'AbsTol', test.tol);    
                    test.verifyEqual( s/2 + 4, y, 'AbsTol', test.tol); 
                case NdgCellType.Quad 
                    test.verifyEqual( r*3 + 2*s + 2, x, 'AbsTol', test.tol);    
                    test.verifyEqual( s/2 + 1, y, 'AbsTol', test.tol); 
            end
        end
        
%         function test_Jacobian(test)
%             switch test.cell.type
%                 case NdgCellType.Line
%                     test.verifyEqual(1/2*ones(size(test.mesh.J)), test.mesh.J, 'AbsTol', test.tol);
%                 case NdgCellType.Tri
%                     test.verifyEqual(1/4*ones(size(test.mesh.J)), test.mesh.J, 'AbsTol', test.tol);                   
%                 case NdgCellType.Quad 
%                     test.verifyEqual(3/2*ones(size(test.mesh.J)), test.mesh.J, 'AbsTol', test.tol); 
%             end
%         end
        
%         function test_transMatrix(test)
%             Dr = test.cell.Dr; rx = test.mesh.rx; ry = test.mesh.sx;
%             Ds = test.cell.Ds; sx = test.mesh.sx; sy = test.mesh.sy;
%             Dt = test.cell.Dt; tx = test.mesh.tx; ty = test.mesh.tx;
%             
%             x = test.mesh.x;
%             y = test.mesh.y;
%             z = test.mesh.z;
%             switch test.cell.type
%                 case NdgCellType.Line
%                     test.verifyEqual(ones(size(test.mesh.rx)), ...
%                         Dr*x.*rx, 'AbsTol', test.tol);
%                 case NdgCellType.Tri
%                     test.verifyEqual(ones(size(test.mesh.rx)), ...
%                         Dr*x.*rx+Ds*x.*sx, 'AbsTol', test.tol);                   
%                     test.verifyEqual(ones(size(test.mesh.rx)), ...
%                         Dr*y.*ry+ Ds*y.*sy, 'AbsTol', test.tol);                     
%                 case NdgCellType.Quad
%                     test.verifyEqual(ones(size(test.mesh.rx)), ...
%                         Dr*x.*rx+Ds*x.*sx, 'AbsTol', test.tol);                   
%                     test.verifyEqual(ones(size(test.mesh.ry)), ...
%                         Dr*y.*ry+Ds*y.*sy, 'AbsTol', test.tol);                   
%             end
%         end
    
    end% methods
end% class

