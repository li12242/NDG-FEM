classdef conv2d_gaussrotate < conv2d
    %CONV2D_GAUSSROTATE Gauss 函数旋转传输算例
    %   Detailed explanation goes here
    
    properties(Constant)
        w = 5*pi/6
        x0 = 0.5
        y0 = 0.5
        r0 = 0.15
        d = 0.25
        theta0 = -pi
    end
    
    methods
        function obj = conv2d_gaussrotate(varargin)
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
                xlim = [0, 1]; ylim = [0, 1]; % 计算域
                Mx = M; My = M; % 单元个数
                zg_bc = ndg_lib.bc_type.ZeroGrad;
                bc_type = [zg_bc, zg_bc, zg_bc, zg_bc];

                input_type = 'uniform';
                input_var = {N, cell_type, xlim, ylim, Mx, My, bc_type};
            end
            
            obj = obj@conv2d(input_type, input_var);
            obj.init; % call initial function
            obj.ftime = 2.4;
        end% func
        
        function init(obj)
            obj.f_Q = obj.ext_func(0);
            % velocity field
            obj.u = obj.w.*(0.5-obj.mesh.y); 
            obj.v = obj.w.*(obj.mesh.x-0.5);
        end% func
        
        
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end
        
        function f_ext = ext_func(obj, time)
            theta = obj.theta0 + obj.w*time; % 旋转角度
            xt = obj.x0 + obj.d*cos(theta);
            yt = obj.y0 + obj.d*sin(theta);
            r2 = sqrt((obj.mesh.x-xt).^2+(obj.mesh.y-yt).^2)./obj.r0;
            ind = ( r2<=1.0);
            
            f_ext = zeros(obj.mesh.cell.Np, obj.mesh.K);
            f_ext(ind) = ( 1+cos(r2(ind)*pi) )./2;
        end
    end
    
end

