classdef swe2d < ndg_lib.phys.phys2d
    %SWE2D 二维浅水方程求解器对象
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 3 % 物理场个数，分别为 h，qx=hu，qy=hv
        gra = 9.80616 % 重力加速度
    end
    
    properties(Abstract, Constant)
        hmin
    end
    
    properties(SetAccess=protected)
        bot % 底坡高程
        bx, by % 底坡高程梯度
        ftime % 计算终止时间
        dt % 计算时间步长
        wetflag % 湿单元逻辑值
        slopelimiter
    end
    
    %% 私有方法
    methods(Access=protected)
        function spe = char_len(obj, f_Q)
            % 计算节点雅克比特征值
            h = f_Q(:,:,1);
            q = sqrt( f_Q(:,:,2).^2 + f_Q(:,:,3).^2 );
            spe = (q./h) + sqrt(obj.gra*h);
            spe(:, ~obj.wetflag) = eps;
        end
        
        function wetdry_detector(obj, f_Q)
            % 判断单元干湿状态
            % 设置所有水深大于阀值的单元类型为湿单元
            obj.wetflag = all( f_Q(:,:,1) > obj.hmin );
            obj.mesh.EToR( ~obj.wetflag ) = ndg_lib.mesh_type.Dry;
            obj.mesh.EToR( obj.wetflag ) = ndg_lib.mesh_type.Normal;
        end
        
        function topo_grad_term(obj)
            % 对于底坡不变的模型，预先计算底坡梯度较少计算量
            obj.bx = obj.mesh.rx.*(obj.mesh.cell.Dr*obj.bot) ...
               + obj.mesh.sx.*(obj.mesh.cell.Ds*obj.bot);
            obj.by = obj.mesh.ry.*(obj.mesh.cell.Dr*obj.bot) ...
               + obj.mesh.sy.*(obj.mesh.cell.Ds*obj.bot);
        end
        
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
        [ dflux ] = surf_term(obj, f_Q ) % 计算单元边界数值通量
        [ E,G ] = flux_term(obj, f_Q ) % 计算通量项
        [ sb ] = topo_sour_term( obj, f_Q ); % 计算底坡源项
        [ dflux ] = hll_surf_term(obj, f_Q ) % 采用 HLL 数值通量计算
        [ dflux ] = lf_surf_term(obj, f_Q ) % 采用 LF 数值通量计算
    end
    
    %% 公共方法
    methods(Abstract)
        [ f_Q ] = init(obj)
    end
    
    methods
        [ obj ] = RK45(obj); % 采用 SSP RK45 格式计算
        [ obj ] = RK45_OBC(obj); % 采用 SSP RK45 格式计算，并考虑开边界条件
        [ obj ] = VB_RK45(obj); % 采用 vertex-based RK 格式计算
        [ obj ] = VB_RK45_OBC(obj); % 采用 VB-RK 格式计算，并考虑开边界条件
        
        function obj = swe2d(mesh)
            obj = obj@ndg_lib.phys.phys2d(mesh);
            %obj.slopelimiter = ndg_utility.limiter.VB.VB_2d(mesh);
            obj.slopelimiter = ndg_utility.limiter.BJ.BJ_2d(mesh);
        end
    end
    
end

