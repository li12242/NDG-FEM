classdef quad_fullquad < ndg_lib.std_cell.quad
    %QUAD_FULLQUADRATURE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nq, Nfq % Gauss-Legrend 体积分与面积分节点个数
        zq, q % 一维 Gauss-Legrend 边界积分节点及权重
        rbq, sbq, wbq % 边界上 Gauss-Legrend 积分节点与权重
        rq, sq, wq  % Gauss-Legrend 积分节点与权重
        Vq, Vbq % 体积分节点与面积分节点投影函数，每列为一个基函数在各点函数值
        Drq, Dsq % 体积分节点处插值函数导数，每列为一个基函数在各点导数
    end
    
    methods
        function obj = quad_fullquad(N)
            obj = obj@ndg_lib.std_cell.quad(N);
            Nq = obj.N+1; % 每个维度积分节点个数
            obj.Nq = Nq^2;
            obj.Nfq = Nq*obj.Nface;
            [obj.zq, obj.q] = Polylib.zwgl(Nq); % 一维 GL 积分节点
            [obj.rq, obj.sq, obj.wq, obj.rbq, obj.sbq, obj.wbq] ...
                = quadrature_point(obj, Nq, obj.zq, obj.q);
            obj.Vq = vandermonde_mat(obj, obj.rq, obj.sq);
            obj.Vbq = vandermonde_mat(obj, obj.rbq, obj.sbq);
            obj.Drq = obj.proj_node2quad(obj.Dr); 
            obj.Dsq = obj.proj_node2quad(obj.Ds);
            %obj.Drq = obj.Drq'; obj.Dsq = obj.Dsq';
        end
        
        function [rq, sq, w, rbq, sbq, wbq] = ...
                quadrature_point(obj, Nq, zq, wq)
            % 单元内 GL 积分节点
            np = numel(zq);
            % 首先按照x坐标循环，然后按照y坐标循环
            rq = zq*ones(1, np);
            sq = ones(np, 1)*zq';

            rq = rq(:); sq = sq(:); w = bsxfun(@times, wq, wq');
            w = w(:);
            % 边界 GL 积分节点
            rbq = zeros(Nq, obj.Nface);
            sbq = zeros(Nq, obj.Nface);
            wbq = zeros(Nq, obj.Nface);
            for f = 1:obj.Nface
                vind = obj.FToV(:, f);
                rv = obj.vr(vind);
                sv = obj.vs(vind);
                rbq(:, f) = 0.5*((1-zq)*rv(1) + (1+zq)*rv(2));
                sbq(:, f) = 0.5*((1-zq)*sv(1) + (1+zq)*sv(2));
                wbq(:, f) = wq;
            end
            rbq = rbq(:); sbq = sbq(:); wbq = wbq(:);
        end
                
        function [ V ] = vandermonde_mat(obj, r, s)
            Ng = numel(r);
            Vg = zeros(Ng, obj.Np);
            for n = 1:obj.Np
                Vg(:, n) = obj.orthogonal_func(obj.N, n, r, s, obj.t);
            end% for
            V = Vg/obj.V;
        end
        
        function [ fq_Q ] = proj_node2quad(obj, f_Q)
            % 根据插值节点系数计算 guass 积分点函数值
            fq_Q = obj.Vq*f_Q;
        end
        
        function [ fbq_Q ] = proj_node2surf_quad(obj, f_Q)
            % 根据插值节点系数计算边界 gauss 积分点函数值
            fbq_Q = obj.Vbq*f_Q;
        end
    end
    
end

