classdef swe2d_Malpasset_dambreak < swe2d
    %swe2d_Malpasset_dambreak Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        hmin = 1e-1
        n = 0.033.^2;
    end
    
    properties(SetAccess=private)
        detector
    end
    
    methods(Access=protected, Hidden)
        function bot = interp_topography(obj) 
            topography_file = ...
                'SWE2D/@swe2d_Malpasset_dambreak/mesh/bathmetry.txt';
            fp = fopen(topography_file);
            fgets(fp);
            data = fscanf(fp, '%e %e %e', [3, inf]);
            fclose(fp);
            interp = scatteredInterpolant(data(1,:)',data(2,:)',...
                data(3,:)','linear');
            bot = interp(obj.mesh.x, obj.mesh.y);
        end
        
        function detector = set_detector(obj)
            % set detector; transformers A B C and gauges 6-14
            xd = [5550, 11900, 13000, 4947.46, 5717.30, 6775.14,...
                7128.20, 8585.30, 9674.97, ...
                10939.15, 11724.37, 12723.70]; 
            yd = [4400, 3250, 2700, ...
                4289.71, 4407.61, 3869.23,...
                3162.00, 3443.08, 3085.89, ...
                3044.78, 2810.41, 2485.08];
            detector = ndg_utility.detector.detector2d(obj.mesh, ...
                xd, yd, 0.5, obj.ftime, obj.Nfield);
        end
    end
    
    methods(Access=protected)
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
        [ sf ] = fric_sour_term( obj, f_Q ) % 计算摩阻源项
    end
    methods
        result(obj)
        [ obj ] = VB_RK45_detect(obj)
        
        function obj = swe2d_Malpasset_dambreak(N, casename, type)
            obj = obj@swe2d(N, casename, type);
            obj.init();
            obj.ftime = 2000;
            
            obj.detector = obj.set_detector();
        end% func
        
        function draw_dam( obj, varargin )
            if nargin == 1
                f_Q = obj.f_Q;
            else
                f_Q = varargin{1};
            end
            obj.mesh.draw( f_Q(:,:,1) ); view([40, 45]); 
            xlim([4200, 5800]); ylim([3600, 4600]);
            drawnow;
        end% func
        
        function init(obj)
            obj.bot = interp_topography(obj);
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            h = zeros(obj.mesh.cell.Np, obj.mesh.K);
            flag = (obj.mesh.EToR == 2);
            h(:, flag) = 100 - obj.bot(:,flag);
%             flag = (obj.mesh.EToR == 1);
%             h(:, flag) = 0 - obj.bot(:, flag);
            h( h< 0 ) = 0;
            f_Q(:, :, 1) = h;
            obj.f_Q = f_Q;
        end% func
    end
end

