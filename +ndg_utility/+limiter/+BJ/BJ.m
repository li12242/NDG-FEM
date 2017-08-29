classdef BJ
    %BJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        cell
        Kv
        VToE
        xc, yc, zc
    end
    
    methods
        function obj = BJ(mesh)
            obj.mesh = mesh;
            obj.cell = mesh.cell;
            [obj.Kv, obj.VToE] = obj.vertex_connect();
            
            [ obj.xc ] = obj.mesh.cell_mean(obj.mesh.x);
            [ obj.yc ] = obj.mesh.cell_mean(obj.mesh.y);
            [ obj.zc ] = obj.mesh.cell_mean(obj.mesh.z);
        end
        
    end
    
    methods(Hidden)
        function [Kv, VToE] =  vertex_connect(obj)
            % Obtain the element index containing each vertex
            Kv = zeros(obj.mesh.Nv, 1);
            for k = 1:obj.mesh.K
                v = obj.mesh.EToV(:, k);
                Kv(v) = Kv(v) + 1;
            end
            maxNe = max(Kv); % The max-number of elements sharing the same vertex
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
    
    methods(Access=protected)
        function [vmax, vmin] = vertex_bound(obj, c_mean)
            % find the maximum and minimum of cell averaged values around 
            % each vertex.
            [ vmax, vmin ] = vertex_extreme(obj.Kv, obj.VToE, c_mean);
        end% func
    end
    
end
