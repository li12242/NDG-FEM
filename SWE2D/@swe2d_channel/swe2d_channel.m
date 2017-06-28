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
    
    properties(SetAccess=private)
        detector
    end
    
    properties
        casename;
        T = 360; % 周期 6 min
        H = 40; % 航道水深
        eta = 0.2; % 潮波振幅
        obc_time_interval = 12; % 开边界数据频率 (s)
    end
    
    methods(Abstract)
        set_bc(obj) % set boundary condition
    end
    
    methods
        draw(obj, f_Q);
        
        function obj = swe2d_channel(N, casename, type)
            [ mesh ] = read_mesh_file(N, casename, type);
            obj = obj@swe2d(mesh);
            obj.casename = casename;
            obj.ftime = 7200;
            
            % set detector
            xd = [0, 5e3, 10e3, 20e3];
            yd = [0, 0, 0, 0];
            obj.detector = ndg_utility.detector.detector2d(mesh, ...
                xd, yd, obj.T/24, obj.ftime, obj.Nfield);
            
            obj.EToB = get_bc_id(obj); % 读取开边界类型
            obj.obc_vert = get_obc_vert(casename); % 读取开边界顶点编号
            
            obj.set_bc(); % 设置开边界条件
        end% func
        
        function init(obj)
            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.f_extQ = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.bot = zeros(obj.mesh.cell.Np, obj.mesh.K);
            % 设定初始条件
            w = 2*pi/obj.T;
            c = sqrt(obj.gra*obj.H);
            k = w/c;
            
            h = obj.eta*cos(k.*obj.mesh.x)+obj.H;
            u = obj.eta*sqrt(obj.gra/obj.H)*cos(k.*obj.mesh.x);

            obj.f_Q = zeros(obj.mesh.cell.Np, obj.mesh.K, obj.Nfield);
            obj.f_Q(:, :, 1) = h;
            obj.f_Q(:, :, 2) = h.*u;
        end
    end
    
    methods(Access=protected)
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
end

