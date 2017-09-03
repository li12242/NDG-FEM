classdef swe2d_dbd < swe2d
    %SWE2D_DBD The dam break problem with dry condition at the downstream
    %   Detailed explanation goes here
    
    properties(Hidden)
        dam_pos = 500
        h0 = 10
    end
    
    properties(Constant)
        hmin = 1e-2
    end
    
    methods
        function obj = swe2d_dbd(varargin)
            [ input ] = set_case_parameter( varargin{:} );
            obj = obj@swe2d(input{:});
            obj.ftime = 20;
            obj.init();
        end% func
        
        function init(obj)
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            xc  = obj.mesh.cell_mean(obj.mesh.x);
            f_Q(:, xc<obj.dam_pos, 1) = obj.h0;
            obj.f_Q = f_Q;
            
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K); % bottom
        end% func
        
        function result(obj)
            % set detector
            np = 100; delta = 1e-3;
            xd = linspace(delta, 1000 -delta, np);
            yd = zeros( size(xd) );
            detector = ndg_utility.detector.detector2d(...
                obj.mesh, xd, yd, obj.ftime/2, obj.ftime, obj.Nfield);
            detector.collect(obj.f_Q, 0);
            % 获取节点上计算解
            h_num = detector.dQ(:, 1, 1); hu_num = detector.dQ(:, 1, 2);
            % 获取节点处精确解
            [h_ext, u_ext] = obj.ext1d_func(xd, obj.ftime);
            figure('color', 'w'); 
            subplot(2,1,1); hold on;
            plot(xd, h_num, 'bo-', 'MarkerSize', 6);
            plot(xd, h_ext, 'r-', 'LineWidth', 2.5);
            box on; grid on;
            
            subplot(2,1,2); hold on;
            plot(xd, hu_num, 'bo-', 'MarkerSize', 6);
            plot(xd, h_ext.*u_ext, 'r-', 'LineWidth', 2.5);
            box on; grid on;
        end
    end
    
    methods(Access=protected)
        function f_ext = ext_func(obj, time)
            f_ext = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            [h, u] = ext1d_func(obj, obj.mesh.x, time);
            
            f_ext(:, :, 1) = h;
            f_ext(:, :, 2) = h.*u;
        end% func
    end
    
    methods(Hidden, Access=private)
        function [h, u] = ext1d_func(obj, x, time)
            h = zeros(size(x));
            u = zeros(size(x));
            % 一维溃坝算例精确解
            temp = (x - obj.dam_pos)/time;
            c0 = sqrt(obj.gra .* obj.h0);
            % left part
            ind = (temp < -c0 );
            h(ind) = obj.h0; u(ind) = 0;
            % middle part
            ind = (temp > -c0 ) & ( temp < 2*c0 );
            h(ind) = (2*c0 - temp(ind) ).^2/obj.gra/obj.gra; 
            u(ind) = 2/3*(c0 + temp(ind) );
        end
    end
    
end

function [ input ] = set_case_parameter( varargin )
% set the range of the computation domain and the boundary conditions
if( isa(varargin{2}, 'char') )
    N = varargin{1};
    casename = varargin{2};
    type = varargin{3};
    input = {N, casename, type};
elseif( isa(varargin{2}, 'double') )
    N = varargin{1};
    M = varargin{2};
    type = varargin{3};
    xlim = [0, 1000]; % computation domain
    ylim = [-10, 10]; % computation domain
    Mx = M; My = 1;
    zg_bc = ndg_lib.bc_type.ZeroGrad;
    bc_type = [zg_bc, zg_bc, zg_bc, zg_bc];
    input = {N, xlim, ylim, Mx, My, type, bc_type};
end% switch

end% func

