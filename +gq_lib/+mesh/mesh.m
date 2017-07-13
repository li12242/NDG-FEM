classdef mesh < ndg_lib.mesh.mesh
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=protected)
        invM % inverse mass matrix of each element
        Jq, % 插值节点处雅克比行列式
        rxwq, rywq, rzwq
        sxwq, sywq, szwq
        txwq, tywq, tzwq
        
        eidMq, eidPq
        eidtypeq,
        nxq, nyq, nzq, Jsq % 面积分系数
    end
    
    methods
        
        function obj = mesh(cell, varargin)
            obj = obj@ndg_lib.mesh.mesh(cell, varargin{:});
            % 计算每个面质量矩阵
            obj.Jq = obj.cell.proj_node2quad(obj.J);
            obj.invM = mass_matrix(obj);
            % 计算积分节点处面积分变换系数与积分权重之积
            obj.rxwq = obj.vol_scal(obj.rx);
            obj.rywq = obj.vol_scal(obj.ry);
            obj.rzwq = obj.vol_scal(obj.rz);
            obj.sxwq = obj.vol_scal(obj.sx);
            obj.sywq = obj.vol_scal(obj.sy);
            obj.szwq = obj.vol_scal(obj.sz);
            obj.txwq = obj.vol_scal(obj.tx);
            obj.tywq = obj.vol_scal(obj.ty);
            obj.tzwq = obj.vol_scal(obj.tz);
            % 边界积分节点处法向量与雅克比系数
            [obj.nxq, obj.nyq, obj.nzq, obj.Jsq] ...
                = face_scal(obj, obj.vx, obj.vy, obj.EToV);
            [obj.eidMq, obj.eidPq, obj.eidtypeq] = connect_quad_node(obj);
        end
    end
    
    methods(Abstract)
        [nxq, nyq, nzq, jsq] = face_scal(obj, vx, vy, EToV);
        [eidMq, eidPq, eidtype] = connect_quad_node(obj);
    end
    
    methods(Hidden, Access=protected)
        function [rxwq] = vol_scal(obj, rx)
            rxwq = bsxfun(@times, obj.cell.wq, ...
                obj.cell.proj_node2quad(rx).*obj.Jq);
        end% func
        
        function invM = mass_matrix(obj)
            invM = zeros(obj.cell.Np, obj.cell.Np, obj.K);
            for k = 1:obj.K % 计算单元质量矩阵
                JVq = bsxfun(@times, obj.cell.wq.*obj.Jq(:, k), ...
                    obj.cell.Vq);
                mass_mat = obj.cell.Vq'*JVq;
                invM(:, :, k) = inv(mass_mat);
            end
        end% func
    end
end

