classdef limiter_vert
    %LIMITER_VERT 基于节点值的斜率限制器。
    %   Detailed explanation goes here
    
    properties
        mesh % 网格对象
        cell % 单元对象
        Kv  % 每个顶点所在单元个数
        VToE % 每个顶点所在单元序号
    end
    
    methods(Hidden)
        function [Kv, VToE] =  vertex_connect(obj)
            % 统计顶点与单元链接关系
            Kv = zeros(obj.mesh.Nv, 1);
            for k = 1:obj.mesh.K
                v = obj.mesh.EToV(:, k);
                Kv(v) = Kv(v) + 1;
            end
            maxNe = max(Kv); % 所有顶点中最多单元数
            VToE = zeros(maxNe, obj.mesh.Nv);
            Kv = zeros(obj.mesh.Nv, 1);
            for k = 1:obj.mesh.K
                v = obj.mesh.EToV(:, k);
                ind = Kv(v)+1 + (v-1)*maxNe;
                VToE(ind) = k;
                Kv(v) = Kv(v) + 1;
            end
        end
    end
    
    methods
        function obj = limiter_vert(mesh, cell)
            obj.mesh = mesh;
            obj.cell = cell;
            [obj.Kv, obj.VToE] = obj.vertex_connect();
        end
    end
    
end

