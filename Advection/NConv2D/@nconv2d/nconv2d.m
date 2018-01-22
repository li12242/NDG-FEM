classdef nconv2d < ndg_lib.phys.phys2d
    %NCONV2D two dimensional nonlinear convection problem.
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 1  % 变量个数
    end
    
    properties
        ftime   % 计算终止时间
        dt      % 计算时间步长
    end
    
    %% 私有函数
    methods(Access=protected) % private 
        [ E, G ] = flux_term( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
        
        function [ spe ] = char_len( obj, f_Q )
            spe = max( max( abs( f_Q ) ) );
        end% func
    end
    
    %% 公共方法
    methods
        function obj = nconv2d(varargin)
            switch varargin{1}
                case 'file' % 读取文件
                    v2 = varargin{2};
                    
                    N = v2{1};
                    cell_type = v2{2};
                    casename = v2{3};
                    [ mesh ] = read_mesh_file(N, cell_type, casename);
                case 'mesh'
                    v2 = varargin{2};
                    mesh = v2{1};
                case 'uniform'
                    v2 = varargin{2};
                    
                    N = v2{1};
                    cell_type = v2{2};
                    xlim = v2{3};
                    ylim = v2{4};
                    Mx = v2{5};
                    My = v2{6};
                    bc_type = v2{7};
                    [ mesh ] = uniform_mesh(N, cell_type, ...
                        xlim, ylim, Mx, My, bc_type);
            end
            
            obj = obj@ndg_lib.phys.phys2d(mesh);
        end% func
        
        f_Q = RK45(obj) % Runge-Kutta 4th order 5 stages
        refine_mesh(obj, multi_ratio) % refined mesh
    end
    
end

