classdef swe1d_trihump < swe1d
    %SWE1D_TRIHUMO Summary of this class goes here
    %   Detailed explanation goes here
    % Reference:
    %   1. Kesserwani & Liang, 2010;
    %   2. Mao et al., 2016
    
    properties(Constant)
        h0 = 0.75
        n = (1.25e-2)^2 % Manning 系数平方
        hmin = 2e-3
    end
    
    properties
        M = 1e-10
        detector
    end
    
    methods(Access=protected)
        [ sf ] = fric_sour_term( obj, f_Q ) % 摩阻源项
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
    
    methods
        function obj = swe1d_trihump(N, K)
            mesh = uniform_mesh(N, K);
            obj = obj@swe1d(mesh);
            
            obj.ftime = 93;
            obj.cfl = 0.2;
            obj.detector = ndg_utility.detector.detector1d(mesh, ...
                15.5+[2, 4, 8, 10, 11, 13, 20], 0.1, obj.ftime, obj.Nfield);
            obj.init();
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
            % detector init
            obj.detector.init();
        end
        
        function result(obj)
            folder = '@swe1d_trihump/private/';
            name = {'p1.txt', 'p2.txt', 'p3.txt', ...
                'p4.txt', 'p5.txt', 'p6.txt', 'p7.txt'};
            for i = 1:numel(name)
                filename = fullfile(folder, name{i});
                figure('color', 'w');
                [t, h] = read_measured_data(filename);
                plot(t, h, 'ro'); hold on;
                p_h = obj.detector.draw(i, 1);
                xlim([0, 90]); ylim([0, .75]);
                grid on; box on;
                xlabel('$time\,(s)$', 'Interpreter', 'Latex', 'FontSize', 16);
                ylabel('$h\,(m)$', 'Interpreter', 'Latex', 'FontSize', 16);
                legend({'Measured', 'NDG'}, 'box', 'off', 'FontSize', 16, ...
                    'Interpreter', 'Latex')
            end
        end
    end
    
end

