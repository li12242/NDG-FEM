classdef conv2d_adv_gq < conv2d
    %CONV2D_ADV_GQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        x0 = -0.5; 
        y0 = -0.5;
        u0 = 0.5; 
        v0 = 0.5;
    end
    
    properties
        uq, vq % flow variables at quadrature points
    end
    
    methods
        function obj = conv2d_adv_gq(varargin)
            if( isa(varargin{2}, 'char') )
                N = varargin{1};
                casename = varargin{2};
                cell_type = varargin{3};
                input_type = 'file';
                input_var = {N, cell_type, casename};
            elseif( isa(varargin{2}, 'double') )
                N = varargin{1}; % 单元阶数
                M = varargin{2}; % 单元个数
                cell_type = varargin{3}; % 单元类型
                xlim = [-1, 1]; ylim = [-1, 1]; % 计算域
                Mx = M; My = M; % 单元个数
                zg_bc = ndg_lib.bc_type.ZeroGrad;
                bc_type = [zg_bc, zg_bc, zg_bc, zg_bc];

                input_type = 'uniform';
                input_var = {N, cell_type, xlim, ylim, Mx, My, bc_type};
            end
            obj = obj@conv2d(input_type, input_var);
            obj.ftime = 2;
            obj.init();
        end% func
        
        function reset_gq_mesh(obj)
            N = obj.mesh.cell.N;
            mesh = obj.mesh;
            switch obj.mesh.cell.type
                case ndg_lib.std_cell_type.Tri
                    stdcell = gq_lib.std_cell.tri(N);
                    obj.mesh = gq_lib.mesh.tri_mesh(stdcell, ...
                        'variable', ...
                        {mesh.Nv, mesh.vx, mesh.vy, ...
                        mesh.K, mesh.EToV, mesh.EToR, mesh.EToBS});
                case ndg_lib.std_cell_type.Quad
                    stdcell = gq_lib.std_cell.quad(N);
                    obj.mesh = gq_lib.mesh.quad_mesh(stdcell, ...
                        'variable', ...
                        {mesh.Nv, mesh.vx, mesh.vy, ...
                        mesh.K, mesh.EToV, mesh.EToR, mesh.EToBS});
            end            
        end% func
        
        function init(obj)
            reset_gq_mesh(obj)
            obj.u = obj.u0*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.v = obj.v0*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.uq = obj.mesh.cell.project_node2quad(obj.u);
            obj.vq = obj.mesh.cell.project_node2quad(obj.u);
            obj.f_Q = obj.ext_func(0);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K);
        end% func
        
    end
    
    %% 私有函数
    methods(Access=protected) % private 
        [ E, G ] = flux_term_quad( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term_quad( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
        
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u*time;
            yc = obj.y0 + obj.v*time;
            
            sigma = 125*1e3/(33*33);
            t = -( (obj.mesh.x-xc).^2+(obj.mesh.y-yc).^2 )*sigma;
            f_ext = exp(t);
        end% func
        
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end% func
    end
    
end

