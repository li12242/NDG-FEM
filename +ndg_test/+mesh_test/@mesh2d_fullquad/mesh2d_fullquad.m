classdef mesh2d_fullquad < ndg_lib.mesh.mesh2d
    %MESH2D_FULLQUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        invM % inverse mass matrix of each element
        Sx, Sy % 体积分刚度矩阵 $S_{x,ij} = J(\xi_j) d\varphi_i/dx|_{\xi_j}$
        Mes % 面积分系数矩阵 Js \cdot \varphi_i(\xi_j)
        Jq, % 插值节点处雅克比行列式
        rxq, ryq, rzq
        sxq, syq, szq
        txq, tyq, tzq
    end
    
    methods
        function obj = mesh2d_fullquad(cell, varargin)
            obj = obj@ndg_lib.mesh.mesh2d(cell, varargin{:});
            % 计算每个面质量矩阵
            obj.Jq = obj.cell.map2vol_quad_point(obj.J);
            obj.invM = mass_matrix(obj);
        end
    
        function invM = mass_matrix(obj)
            invM = zeros(obj.cell.Np, obj.cell.Np, obj.K);
            for k = 1:obj.K % 计算单元质量矩阵
                JVq = bsxfun(@times, obj.cell.wq.*obj.Jq(:, k), ...
                    obj.cell.Vq);
                mass_mat = obj.cell.Vq'*JVq;
                invM(:, :, k) = inv(mass_mat);
            end
        end
        
        function [ Sx, Sy ] = stiff_matrix(obj)
            
        end
    end
end

