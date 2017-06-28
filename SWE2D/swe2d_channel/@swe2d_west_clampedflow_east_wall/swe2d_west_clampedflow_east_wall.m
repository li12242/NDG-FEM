classdef swe2d_west_clampedflow_east_wall < swe2d_channel
    %SWE2D_WEST_CLAMPED_EAST_RADIATION Summary of this class goes here
    %   Detailed explanation goes here
        
    properties
        flux = 20; % m2/s
    end
    methods
        function obj = swe2d_west_clampedflow_east_wall(N)
            casename = 'SWE2D/@swe2d_channel/mesh/quad';
            type = ndg_lib.std_cell_type.Quad;
            obj = obj@swe2d_channel(N, casename, type);
        end
        
        function set_bc(obj)
            % 设定开边界条件 EToBS
            EToBS = obj.mesh.EToBS;
            westid = (obj.EToB == 2); 
            EToBS(westid) = ndg_lib.bc_type.ClampedVel; % 西侧边界条件
            eastid = (obj.EToB == 3);
            EToBS(eastid) = ndg_lib.bc_type.SlipWall; % 东侧边界条件
            % 生成新网格对象
            obj.mesh = ndg_lib.mesh.mesh2d(obj.mesh.cell, ...
                obj.mesh.Nv, obj.mesh.vx, obj.mesh.vy, ...
                obj.mesh.K, obj.mesh.EToV, obj.mesh.EToR, EToBS);

            % 设定边界信息
            Nv = numel(obj.obc_vert); % 顶点个数
            Nt = 1;
            Nfield = 3; % 开边界物理场个数
            time = 0; % 开边界数据时间
            f_Q = zeros(Nv, Nfield, Nt);
            f_Q(:, 1, :) = obj.H;
            f_Q(:, 2, :) = obj.flux;
            
            % 生成开边界文件
            obj.obc_file = ndg_lib.phys.obc_file();
            obj.obc_file.make_obc_file([obj.casename, '.nc'], ...
                time, obj.obc_vert, f_Q);
            obj.obc_file.set_file([obj.casename, '.nc']);

            % 设置初始条件
            obj.init();
        end
        
        function init(obj)
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
            obj.f_Q(:,:,1) = obj.H;
        end
    end
    
end

