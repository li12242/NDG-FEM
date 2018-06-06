classdef conv2d_solidrotate < conv2d
    %CONV2D_SOLIDROTATE Summary of this class goes here
    %   Detailed explanation goes here
    % Reference:
    %   1. Kuzmin, 2014;
    
    properties(Constant)
        limitertype = 'WENO'
    end
    properties
        slopelimiter
        M = 1e-3;
    end
    
    methods
        function [ spe ] = character_len(obj, f_Q)
            vel = sqrt( obj.u.^2 + obj.v.^2 );
            spe = vel;
        end
        
        function obj = conv2d_solidrotate(varargin)
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
            obj.init(); % call initial function
            obj.ftime = 2.4;
            obj.slopelimiter = ndg_utility.limiter.TVB.TVB_tri(obj.mesh);
        end
        
        function init(obj)
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K);
            r0  = 0.15;
            % The slotted cylinder
            x0  = 0.5; y0  = 0.75;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0) & ...
                ( (abs(obj.mesh.x - x0)>= 0.025) | (obj.mesh.y>=0.85) );
            obj.f_Q(ind) = 1;
            % The cone
            x0 = 0.5; y0 = 0.25;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0);
            obj.f_Q(ind) = 1 - r2(ind);
            % The hump
            x0 = 0.25; y0 = 0.5;
            r2  = sqrt((obj.mesh.x-x0).^2+(obj.mesh.y-y0).^2)./r0;
            ind = ( r2<=1.0);
            obj.f_Q(ind) = (1 + cos(r2(ind)*pi))*0.25;
            % velocity field
            w = 5*pi/6;
            obj.u = w.*(0.5-obj.mesh.y); 
            obj.v = w.*(obj.mesh.x-0.5);
        end
    end
    
end

