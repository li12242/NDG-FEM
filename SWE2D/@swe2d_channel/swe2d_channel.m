classdef swe2d_channel < swe2d
    %SWE2D_CHANNEL 模拟潮波在航道内流动。
    %   测试不同的开边界条件，分别为：
    %       1. 西侧0梯度，东侧水位指定边界条件（未给出外部值的变量皆为0梯度边界条件）；
    %       2. 西侧0梯度，东侧 clamped 边界条件；
    %       3. 西侧指定流量，东侧可滑移固壁；
    
    properties(Constant)
        hmin = 1e-2
    end
    
    properties(Hidden)
        EToB; % 边界序号：0-内部、1-南侧与北侧、2-西侧、3-东侧；
        obc_vert; % 开边界顶点序号
    end
    
    properties
        casename;
        T = 360; % 周期 6 min
        H = 40; % 航道水深
        eta = 0.2; % 潮波振幅
        obc_time_interval = 12; % 开边界数据频率 (s)
    end
    
    methods
        draw(obj, f_Q);
        set_west_wave_depth_obc1(obj, casename);
        set_west_wave_all_obc2(obj, casename);
        set_west_wave_flow_east_wall_obc3(obj, casename)
        set_east_spg_obc4(obj, casename);
        
        function obj = swe2d_channel(N, casename, type)
            [ mesh ] = read_mesh_file(N, casename, type);
            obj = obj@swe2d(mesh);
            obj.casename = casename;
            obj.ftime = 7200;
            obj.set_west_wave_depth_obc1(); % default obc 
        end% func
        
        function init(obj)
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
        end% func
    end
    
    methods(Access=protected)
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
end

