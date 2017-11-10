classdef StdSubCellTest < matlab.unittest.TestCase
    %STDSUBCELLTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> test cell types
        type = {...
%             NdgCellType.Line, ...
            NdgCellType.Tri, ...
%             NdgCellType.Quad
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
            switch type
                case NdgCellType.Tri
                    test.cell = StdSubTri( order );
            end
        end% func
    end
    
    methods(Test, ParameterCombination = 'sequential')
        function TestEdgeNormalVector( obj )
            for n = 1:obj.cell.NFV
            end
        end
        
        function TestProjectMatrix( obj )
        end
    end% methods
    
end

