classdef BJ
    %BJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        cell
    end
    
    methods
        function obj = BJ(mesh)
            obj.mesh = mesh;
            obj.cell = mesh.cell;
        end
    end
    
end

