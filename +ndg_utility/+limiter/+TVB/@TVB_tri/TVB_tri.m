classdef TVB_tri < ndg_utility.limiter.TVB.TVB
    %TVB_TRI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xf, yf % 单元边界中心点坐标
    end
    
    methods
        function obj = TVB_tri(mesh)
            obj = obj@ndg_utility.limiter.TVB.TVB(mesh, mesh.cell);
            
            if (mesh.cell.type ~= ndg_lib.std_cell_type.Tri)
                error('The input cell type is incorrect!');
            end
            obj.xf = mesh.face_mean(mesh.x);
            obj.yf = mesh.face_mean(mesh.y);
            obj.xc = mesh.cell_mean(mesh.x);
            obj.yc = mesh.cell_mean(mesh.y);
        end
        
        function f_limQ = limit(obj, f_Q, M)
            f_cQ = obj.mesh.cell_mean(f_Q);
            f_fQ = obj.mesh.face_mean(f_Q);
            f_limQ = TVB(f_Q, obj.mesh.x, obj.mesh.y, ...
                f_cQ, obj.xc, obj.yc, f_fQ, obj.xf, obj.yf, ...
                obj.mesh.EToE, M);
        end
    end
    
end

