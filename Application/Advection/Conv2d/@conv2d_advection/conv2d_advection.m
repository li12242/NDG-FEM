classdef conv2d_advection < conv2d
    %CONV2D_ADVECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        x0 = -0.5; 
        y0 = -0.5;
        u0 = 0.5;
        v0 = 0.5;
    end
    
    %% private methods
    methods(Access=protected)
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end% func
    end
    
    methods(Access=protected)
        function f_ext = ext_func(obj, time)
            % 模拟精确解
            xc = obj.x0 + obj.u*time;
            yc = obj.y0 + obj.v*time;
            
            sigma = 125*1e3/(33*33);
            t = -( (obj.mesh.x-xc).^2+(obj.mesh.y-yc).^2 )*sigma;
            f_ext = exp(t);
        end% func
    end
    
    %% public methods
    methods
        obj = RK45_section(obj);

        function obj = conv2d_advection(varargin)
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
                cl_bc = ndg_lib.bc_type.Clamped;
                bc_type = [cl_bc, cl_bc, cl_bc, cl_bc];

                input_type = 'uniform';
                input_var = {N, cell_type, xlim, ylim, Mx, My, bc_type};
            end
            obj = obj@conv2d(input_type, input_var);
            obj.ftime = 2;
            obj.init();
        end% func
                
        function init(obj)
            one = ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.u = obj.u0*one;
            obj.v = obj.v0*one;
            obj.f_Q = obj.ext_func(0);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K);
        end% func
    end
    
end

