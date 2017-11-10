classdef conv1d_advection < conv1d
    %CONV1D_ADVECTION Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        L = 2*pi/1.5;
        c = 0.5;
        ob_frec = 200; % number of obc values
    end
    properties(SetAccess=protected)
        slopelimiter
    end
    
    methods
        function obj = conv1d_advection(varargin)
            
            switch nargin
                case 1
                    mesh = varargin{1};
                case 2
                    N = varargin{1};
                    K = varargin{2};
                    mesh = uniform_mesh(N, K); % 调用私有函数构造网格
                otherwise
                    error('The number of input variable is incorrect!');
            end
            
            obj = obj@conv1d(mesh);
            obj.init; % call initial function
            obj.cfl = 0.3;
            obj.ftime = 1.5;
            obj.obc_file = obj.set_obc();
            obj.slopelimiter = ndg_utility.limiter.BJ(mesh);
        end% func
        
        function init(obj) 
            % set the initial value
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K);
            obj.u = obj.c*ones(obj.mesh.cell.Np, obj.mesh.K);
            obj.f_Q = obj.ext_func(0);
        end% func
                
        function [ spe ] = character_len(obj, f_Q)
            vel = abs( obj.u );
            spe = vel;
        end
        
        function update_ext(obj, time)
            obj.f_extQ = ext_func(obj, time);
        end
        
        function obc_file = set_obc(obj)
            vert = [1, obj.mesh.Nv];
            vx = obj.mesh.vx(vert);
            time = linspace(0, obj.ftime, obj.ob_frec);
            Nv = 2; Nt = obj.ob_frec; Nfield = 1;
            f_extQ = zeros(Nv, Nfield, Nt);
            for t = 1:Nt
                tloc = time(t);
                f_extQ(:, 1, t) = sin( (vx - obj.c*tloc)*obj.L );
            end
            
            obc_file = ndg_lib.phys.obc_file();
            filename = './Conv1D/@conv1d_advection/conv1d_advection.nc';
            obc_file.make_obc_file(filename, time, vert, f_extQ);
            obc_file.set_file(filename);
        end
        
        function [ f_ext ] = ext_func(obj, time)
            f_ext = sin( (obj.mesh.x - obj.c*time)*obj.L );
        end
    end
    
end

