classdef conv2d_gaussquad < conv2d
    %CONV2D_GAUSSQUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        x0 = -0.5; y0 = -0.5;
        u0 = 1/2; v0 = 1/2;
        miu = 0;
    end
    
    properties
        uq, vq
    end
    
    methods
        function obj = conv2d_gaussquad(varargin)
            if( isa(varargin{2}, 'char') )
                N = varargin{1};
                casename = varargin{2};
                type = varargin{3};
                [ mesh ] = read_mesh_file(N, casename, type);
            elseif( isa(varargin{2}, 'double') )
                N = varargin{1};
                M = varargin{2};
                type = varargin{3};
                mesh = uniform_mesh(N, M, type);
            end
            obj = obj@conv2d(mesh);
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
            obj.uq = obj.mesh.cell.proj_node2quad(obj.u);
            obj.vq = obj.mesh.cell.proj_node2quad(obj.u);
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
    
    %% Ë½ÓÐº¯Êý
    methods(Access=protected) % private 
        [ E, G ] = flux_term_quad( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term_quad( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
    end
    
end

