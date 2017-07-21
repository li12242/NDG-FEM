classdef conv2d_diffusion < conv2d
    %CONV2D_DIFFUSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        x0 = -0.5; 
        y0 = -0.5;
        u0 = 0.5;
        v0 = 0.5;
        miu = 0;
    end
    
    methods
        obj = RK45_section(obj);
        
        function obj = conv2d_diffusion(varargin)
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
        end
        
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end
        
        function init(obj)
            obj.u = obj.u0*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.v = obj.v0*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.f_Q = obj.ext_func(0);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K);
            
        end
        
        function f_ext = ext_func(obj, time)
            xc = obj.x0 + obj.u*time;
            yc = obj.y0 + obj.v*time;
            if obj.miu > 0
                t = -(obj.mesh.x-xc).^2/obj.miu ...
                    -(obj.mesh.y-yc).^2/obj.miu;
            else
                sigma = 125*1e3/(33*33);
                t = -( (obj.mesh.x-xc).^2+(obj.mesh.y-yc).^2 )*sigma;
            end
            f_ext = exp(t);
        end
    end
    
end

