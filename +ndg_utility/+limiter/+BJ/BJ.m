classdef BJ
    %BJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        cell
    end
    
    methods
        function obj = BJ(mesh, cell)
            obj.mesh = mesh;
            obj.cell = cell;
        end
    end
    
end

