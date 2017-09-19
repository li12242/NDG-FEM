%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef UMeshTest < matlab.unittest.TestCase
    %UMESHTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(MethodSetupParameter)
        types = {...
            StdCellType.Line, ...
            StdCellType.Tri, ...
            StdCellType.Quad, ...
            }
        order = {1, 2, 3}
    end
    
    properties
        mesh
    end
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_std_cell(test, type, order)
            test.mesh = setTestMesh(order, type);
        end% func
    end
    
    methods
    end
    
end

function setTestMesh(order, type)
switch type
    case StdCellType.Line
    case StdCellType.Tri
    case StdCellType.Quad
end% switch
end% func

function setUnionLineMesh()
end

function setUnionTriMesh()
end

function setUnionQuadMesh()
end

function setMixMesh2d()
end

