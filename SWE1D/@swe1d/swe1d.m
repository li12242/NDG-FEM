classdef swe1d < ndg_lib.phys.phys1d
    %SWE1D 一维浅水方程求解器
    %   一维守恒浅水方程，以水深 h，流量 q = hu 为自变量，方程包括连续方程与动量方程
    %       dh/dt + dq/dx = 0
    %       dq/dt + d(uh^2 + 0.5gh^2)/dx = Sb + Sf
    %   其中 Sb 与 Sf 分别为底坡源项与摩阻源项，其中 Sb = -ghdz/dx，摩阻源项
    %   有 Manning 公式与线性公式两种形式。
    %
    properties(Constant)
        Nfield = 2  % 物理量场个数
        gra = 9.81  % 重力加速度
    end
    properties(Abstract, Constant)
        hmin % 最小水深阀值
    end
    
    properties(Abstract)
        M   % TVB 限制器系数
    end

    properties(SetAccess=protected)
        bot     % 底坡高程
        cfl     % CFL 数
        ftime   % 计算终止时间
        dt      % 计算时间步长
        wetflag % 湿单元逻辑数组
        slopelimiter % 斜率限制器
    end

    %% 虚函数
    methods(Abstract)
        init(obj) % 初始化函数
    end
    
    %% 私有函数
    methods(Access=protected)
        [ E ] = flux_term( obj, f_Q ) % 计算体积积分通量项 F
        [ dflux ] = surf_term( obj, f_Q )
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
        [ S ] = topo_sour_term( obj, f_Q ) % 底坡源项 
        [ f_Q ] = positive_preserve( obj, f_Q )
        
        function wetdry_detector(obj, f_Q)
            %hm = obj.mesh.cell_mean(f_Q(:,:,1)); % 计算单元平均水深
            obj.wetflag = all(f_Q(:,:,1) > obj.hmin);
            %obj.wetflag = all( f_Q(:,:,1) > obj.hmin );
            % 设置平均水深小于阀值的单元类型为干单元
            obj.mesh.EToR( ~obj.wetflag ) = ndg_lib.mesh_type.Dry;
            obj.mesh.EToR( obj.wetflag ) = ndg_lib.mesh_type.Normal;
        end
        
        function dt = time_interval( obj, f_Q )
            h = f_Q(:,:,1);
            q = f_Q(:,:,2);
            u = abs(q./h) + sqrt(obj.gra*h);
            s = bsxfun(@times, obj.mesh.vol/obj.mesh.cell.N, 1./u);
            dt = obj.cfl*min( min( s(:, obj.wetflag) ) );
        end
    end

    %% 公共函数
    methods
        [ dflux ] = hll_surf_term( obj, f_Q ) % 计算边界积分通量差值 (Fn - Fn*)
        [ dflux ] = roe_surf_term( obj, f_Q )
        
        function obj = swe1d(mesh)
            obj = obj@ndg_lib.phys.phys1d(mesh);
            obj.slopelimiter = ndg_utility.limiter.TVB(mesh);
        end
    end

end
