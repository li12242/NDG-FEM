classdef phys < matlab.mixin.SetGet
    %PHYS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        Nfield  % 变量个数
    end
    
    properties(Hidden=true)
        draw_h  % 立体绘制图像句柄
    end
    
    properties(SetAccess=protected)
        mesh        % 网格对象
        f_extQ      % 外部值
        obc_file    % 开边界文件
        out_file    % 结果文件
    end
    properties
        f_Q     % 变量
    end
    %% 虚函数
    methods(Abstract)
        draw(obj, field) % 绘制场图
        [ f_Q ] = init(obj) % 初始化
    end
    
    methods(Abstract, Access=protected)
        [ E ] = flux_term( obj, f_Q ) % 计算体积积分通量项 F
        [ dflux ] = surf_term( obj, f_Q ) % 计算边界积分通量差值 (Fn - Fn*)
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
    
    %% 公共方法
    methods
        function obj = phys(mesh)
            obj.mesh = mesh;
            obj.f_Q = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
            obj.f_extQ = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
        end% func
    end
    
    % 非恒定问题时间离散方法
    methods
        f_Q = RK45_solve(obj) % Runge-Kutta 4th order 5 stages
    end
    
    % 文件 I/O
    methods
        function obj = set_out_file(obj, filename, dt)
            % 设置输出文件对象
            obj.out_file = ndg_lib.phys.out_file(filename, ...
                obj.mesh.cell.Np, obj.mesh.cell.K, dt);
        end% func
        
        function obj = set_obc_file(obj, filename)
            % 设置开边界文件对象
            obj.obc_file = ndg_lib.phys.obc_file(filename);
        end% func
        
        function obj = update_ext(obj, stime)
            % 根据开边界文件结果更新外部数据
            vert_extQ = obj.obc_file.get_extQ(obj.Nfield, stime);
            vertlist = obj.obc_file.vert;
            vert_Q = zeros(obj.mesh.Nv, obj.Nfield);
            for fld = 1:obj.Nfield
                vert_Q(vertlist + (fld-1)*obj.mesh.Nv ) = vert_extQ(:, fld);
            end
            obj.f_extQ = obj.mesh.proj_vert2node(vert_Q);
        end% func
        
        function obj = init_from_file(obj, filename)
            % 读取文件数据进行初始化
            fp = fopen(filename);
            Num = fscanf(fp, '%d', 1);
            Nfld = fscanf(fp, '%d', 1); % read number of physical fields
            if ( ( (Num~=obj.mesh.K) && (Num~=obj.mesh.Nv) ) )
                error(['The number of values in file: ', ...
                    num2str(Num), ...
                    ' is neither element number: ', num2str(obj.mesh.K), ...
                    ' nor vertex number: ', num2str(obj.mesh.Nv)]);
            elseif (Nfld~=obj.Nfield)
                error(['The number of physical field in file: ', ...
                    num2str(Nfld), ...
                    ' is different from this phys object: ', ...
                    num2str(obj.Nfield)]);
            end
            fmtStr = ['%d ', repmat('%g ', 1, Nfld)];
            data = fscanf(fp, fmtStr, [Nfld+1, Num]);
            switch Num
                case obj.mesh.K
                    fprintf('\nInit with elemental averaged values.\n\n')
                case obj.mesh.Nv
                    fprintf('\nInit with vertex values.\n\n')
            end
            fclose(fp);
        end
    end
    
end