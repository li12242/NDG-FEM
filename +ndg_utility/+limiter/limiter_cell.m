classdef limiter_cell
    %LIMITER_CELL 基于单元梯度限制的斜率限制器
    %   单元斜率限制器根据周围单元均值对每个单元内自变量的梯度进行抑制，并且根据限制后
    %   的梯度进行重构，获得修正后的单元值。
    
    properties
        phys
        mesh
        cell
    end
    
    methods
        function obj = limiter_cell(mesh, cell)
            obj.mesh = mesh;
            obj.cell = cell;
        end
    end
    
end

