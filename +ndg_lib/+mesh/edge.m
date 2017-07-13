classdef edge
    %EDGE 网格内边界对象
    %   Detailed explanation goes here
    
    properties
        cell % 标准单元
    end
    
    properties(SetAccess=protected)
        Nedge % 边界个数
        Nnode % 边界节点个数
        kM, kP % 相邻两侧单元编号
        fM, fP % 所在相邻两个单元面编号
        ftype@int8 % 面类型
        idM, idP % 边界上两侧节点编号
        fpM, fpP % 边界上两侧节点在面节点中局部编号
        fscal % 雅克比行列式
        fnxM, fnyM, fnzM % 法向量
    end
    
    methods
        function obj = edge(cell, x, y, z, EToV, EToE, EToF, EToBS)
            
            obj.cell = cell;
            [obj.Nedge, obj.Nnode, obj.kM, obj.kP, obj.fM, obj.fP, obj.ftype] ...
                = obj.edge_connect(EToV, EToE, EToF, EToBS);
            
        end
    end
    
    methods(Access=protected)
        function [ Nedge, kM, kP, fM, fP, ftype ] ...
                = edge_connect(obj, EToV, EToE, EToF, EToBS)
            
        end
        
        function [Nnode, iM, iP, fpM, fpP] ...
                = node_connect(obj, x, y, z, Nedge, kM, kP, fM, fP)
        end
        
        
    end
end

