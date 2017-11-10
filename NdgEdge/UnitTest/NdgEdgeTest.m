classdef NdgEdgeTest < matlab.unittest.TestCase

    properties(MethodSetupParameter)
        %> test edge types
        type = { ...
            %'NdgEdge1d', ...
            'NdgEdge2d', ...
            }
        %> test cell orders
        order = {3,4,5}
    end
    
    properties(Constant)
        %> tolerance
        tol = 1e-9;
        %> mixed file name
        mixfile2d = 'mesh2d/MixMesh.msh';
    end
    
    properties
        meshUnion
    end
    
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_test_edge(test, type, order)
            switch type
                case 'NdgEdge1d'
                case 'NdgEdge2d'
                    test.meshUnion = makeGmshFileUMeshUnion2d( order, test.mixfile2d );
            end% switch
        end% func
    end
    
    methods(Test, ParameterCombination = 'sequential')
        function testEdgeConnect( test )
            mesh = test.meshUnion;
            for n1 = 1:numel(mesh)
                Nedge = numel(mesh(n1).edgeUnion);
                for e = 1:Nedge
                    edge = mesh(n1).edgeUnion(e);
                    n2 = edge.FToM;
                    for m = 1:edge.M
                        k1 = edge.FToE(1, m);
                        k2 = edge.FToE(2, m);
                        f1 = edge.FToF(1, m);
                        f2 = edge.FToF(2, m);
                        loc_v1 = mesh(n1).cell.FToV(:, f1);
                        loc_v2 = mesh(n2).cell.FToV(:, f2);
                        v1 = sort( mesh(n1).EToV(loc_v1, k1) );
                        v2 = sort( mesh(n2).EToV(loc_v2, k2) );
                        test.verifyEqual(v1, v2, 'AbsTol', test.tol);
                    end
                end
            end
        end
        
        function testInterpolationMatrix( test )
            mesh = test.meshUnion;
            for n1 = 1:numel(mesh)
                Nedge = numel( mesh(n1).edgeUnion );
                for e = 1:Nedge
                    edge = mesh(n1).edgeUnion(e);
                    n2 = edge.FToM;
                    IntMat = edge.IntM;
                    for m = 1:edge.M
                        k1 = edge.FToE(1, m);
                        k2 = edge.FToE(2, m);
                        p1 = edge.FToN1(:, m);
                        p2 = edge.FToN2(:, m);
                        
                        x1 = mesh(n1).x(p1,k1);
                        x2 = mesh(n2).x(p2,k2);
                        temp = IntMat*x2;
                        test.verifyEqual( x1, temp, 'AbsTol', test.tol);
                    end
                end
            end
        end% func
    end
end

