classdef swe1d_trihump < swe1d
    %SWE1D_TRIHUMO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        h0 = 0.75
        n = (1.25e-2)^2 % Manning 系数平方
        hmin = 2e-3
    end
    
    properties
        M = 1e-10
    end
    
    methods(Access=protected)
        [ sf ] = fric_sour_term( obj, f_Q ) % 摩阻源项
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
    
    methods
        function obj = swe1d_trihump(N, K)
            mesh = uniform_mesh(N, K);
            obj = obj@swe1d(mesh);
            obj.init();
            obj.ftime = 93;
            obj.cfl = 0.2;
        end
        
        function init(obj)
            % initial condition
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            xc = obj.mesh.cell_mean(obj.mesh.x); dam_pos = 15.5;
            f_Q(:, xc< dam_pos, 1) = obj.h0;
            obj.f_Q = f_Q;
            % bottom topography
            bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
            x0 = 28.5; len = 6/2; b0 = 0.4;
            ind = (obj.mesh.x > (x0-len)) & (obj.mesh.x < (x0+len));
            bot(ind) = b0 - abs(obj.mesh.x(ind) - x0)*b0/len;
            obj.bot = bot;
        end
    end
    
end

