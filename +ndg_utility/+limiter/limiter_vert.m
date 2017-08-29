classdef limiter_vert
    %LIMITER_VERT The vertex-based limiter
    %   Detailed explanation goes here
    
    properties
        mesh % mesh object
        cell % cell object
        Kv  % number of elements for each vertex belongs to
        VToE % index of elements contain each vertex
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

