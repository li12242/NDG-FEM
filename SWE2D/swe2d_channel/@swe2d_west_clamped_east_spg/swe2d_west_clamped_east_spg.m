classdef swe2d_west_clamped_east_spg < swe2d_channel
    %SWE2D_WEST_CLAMPED_EAST_SPG Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = swe2d_west_clamped_east_spg(N)
            casename = 'SWE2D/@swe2d_channel/mesh/quad_sponge';
            type = ndg_lib.std_cell_type.Quad;
            obj = obj@swe2d_channel(N, casename, type);
        end
        
        function set_bc(obj)
            % 设定开边界条件 EToBS
            EToBS = obj.mesh.EToBS;
            westid = (obj.EToB == 2); 
            EToBS(westid) = ndg_lib.bc_type.Clamped; % 西侧边界条件
            eastid = (obj.EToB == 3);
            EToBS(eastid) = ndg_lib.bc_type.ZeroGrad; % 东侧边界条件
            % 生成新网格对象
            obj.mesh = ndg_lib.mesh.mesh2d(obj.mesh.cell, ...
                obj.mesh.Nv, obj.mesh.vx, obj.mesh.vy, ...
                obj.mesh.K, obj.mesh.EToV, obj.mesh.EToR, EToBS);
            
            % 添加海绵层
            obj.mesh = obj.mesh.add_sponge(obj.obc_vert);

            % 设定边界信息
            Nt = ceil(obj.ftime/obj.obc_time_interval);
            Nfield = 3; % 开边界物理场个数
            spg_vert = unique(  ... % 海绵层内节点
                obj.mesh.EToV(:, ...
                obj.mesh.EToR == ndg_lib.mesh_type.Sponge) ); 
            vert = unique([spg_vert; obj.obc_vert]); % 全部边界节点
            Nv = numel(vert); % 顶点个数
            vx = obj.mesh.vx(vert);
            time = linspace(0, obj.ftime, Nt); % 开边界数据时间
            f_Q = zeros(Nv, Nfield, Nt);

            % 开边界数据
            w = 2*pi/obj.T;
            c = sqrt(obj.gra*obj.H);
            k = w/c;
            for t = 1:Nt
                tloc = time(t);
                temp = cos(k.*vx - w*tloc);
                h = obj.eta*temp + obj.H;
                u = obj.eta*sqrt(obj.gra/obj.H)*temp;
                f_Q(:, 1, t) = h;
                f_Q(:, 2, t) = h.*u;
            end
            
            for n = 1:Nv
                vid = vert(n);
                if any(spg_vert == vid)
                    f_Q(n, 1, :) = obj.H;
                    f_Q(n, 2, :) = 0;
                end
            end

            % 生成开边界文件
            obj.obc_file = ndg_lib.phys.obc_file();
            obj.obc_file.make_obc_file([obj.casename, '.nc'], ...
                time, vert, f_Q);
            obj.obc_file.set_file([obj.casename, '.nc']);
            
            % 设置初始条件
            obj.init();
        end     
        
    end
    
end

