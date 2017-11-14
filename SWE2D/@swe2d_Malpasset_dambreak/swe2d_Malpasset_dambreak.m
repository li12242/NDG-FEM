classdef swe2d_Malpasset_dambreak < swe2d
    %SWE2D_TSUAMI Summary of this class goes here
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
                 'SWE2D/@swe2d_Malpasset_dambreak/mesh/bathmetry1.txt';
            fp = fopen(topography_file);
            fgets(fp);
            data = fscanf(fp, '%e %e %e', [3, inf]);
            fclose(fp);
            interp = scatteredInterpolant(data(1,:)',data(2,:)',...
                data(3,:)','linear');
            bot = interp(obj.mesh.x, obj.mesh.y);
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
            [ mesh ] = read_mesh_file(N, casename, type);
            obj = obj@swe2d(mesh);
            obj.init();
            obj.ftime = 2000;
%             obj.obc_file = obj.make_obc_file(casename);
            % set detector; transformers A B C ,gauges 6-14, and points 1-17
            xd = [5550, 11900, 13000, 4947.46, 5717.30, 6775.14,...
                7128.20, 8585.3, 9674.97, 10939.15, 11724.37, 12723.70...
                4913.11 5159.75 5790.63 5886.54 6763.05 6929.97 7326.02 7441.01...
                8735.94 8628.6 9761.13 9800 10957 11156.99 11689.05 11626.05 12333.72]; 
            yd = [4400, 3250, 2700, 4289.71, 4407.61, 3869.23,...
                3162.00, 3443.08, 3085.89, 3044.78, 2810.41, 2485.08...
                4244.01 4369.62 4177.76 4503.97 3429.6 3591.87 2948.78 3232.12 3264.61...
                3604.63 3480.36 2414.79 2651.94 3800.72 2592.36 3406.8 2269.74];
            obj.detector = ndg_utility.detector.detector2d(mesh, ...
                xd, yd, 0.5, obj.ftime, obj.Nfield);
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
        
        function obj = update_ext(obj, stime)
            % 根据开边界文件结果更新外部数据,obc_file为double类
            vert_extQ = obj.obc_file.get_extQ(stime);
            vertlist = obj.obc_file.vert;
            if (stime < obj.obc_file.time(end))
                vert_Q = zeros(obj.mesh.Nv, 1);
                vert_Q( vertlist ) = vert_extQ(:, 1); % 获取水位外部值
                obj.f_extQ(:,:,1) = obj.mesh.proj_vert2node(vert_Q) - obj.bot;

                vert_Q = zeros(obj.mesh.Nv, 1);
                vert_Q( vertlist ) = vert_extQ(:, 2); % 获取流量外部值
                obj.f_extQ(:,:,2) = obj.mesh.proj_vert2node(vert_Q);

                vert_Q = zeros(obj.mesh.Nv, 1);
                vert_Q( vertlist ) = vert_extQ(:, 3); % 获取流量外部值
                obj.f_extQ(:,:,3) = obj.mesh.proj_vert2node(vert_Q);
            else
                obj.f_extQ = obj.f_Q;
            end
        end% func
        
        function init(obj)
            obj.bot = interp_topography(obj);
            f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            h = zeros(obj.mesh.cell.Np, obj.mesh.K);
            
            Coor_dam_x = [4701.18,4656.5];Coor_dam_y = [4143.41,4392.1];
            slope = (Coor_dam_y(2)-Coor_dam_y(1))/(Coor_dam_x(2)-Coor_dam_x(1));
            
            flag = obj.mesh.y - Coor_dam_y(1) - slope*(obj.mesh.x - Coor_dam_x(1));
            
            I = find( flag <= 0 );
            
            h(I) = 100 - obj.bot(I);
            
            J = ( obj.mesh.x >= 4394 ).* ( obj.mesh.x <= 5000 ).* ...
                ( obj.mesh.y >= 4550 ).* ( obj.mesh.y <= 5500 ) ;
            
            h( J == 1 ) = 0;
            
%             flag = (obj.mesh.EToR == 2);
%             h(:, flag) = 100 - obj.bot(:,flag);
% %             flag = (obj.mesh.EToR == 1);
% %             h(:, flag) = 0 - obj.bot(:, flag);
            h( h< 0 ) = 0;
% %             flag = find(obj.mesh.EToR == 1);
% %             h_temp = zeros(obj.mesh.cell.Np, obj.mesh.K);
% %             h_temp(:,flag) = obj.bot(:,flag);
% %             Index = find(h_temp<0);
% %             h_temp(Index) = -obj.bot(Index);
% %             h(:,flag) =( h_temp(:,flag)-obj.bot(:,flag) )/2;

            f_Q(:, :, 1) = h;
            obj.f_Q = f_Q;
        end% func
    end
end

