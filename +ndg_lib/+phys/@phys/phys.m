classdef phys < handle
    %PHYS physical field
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        Nfield  % No. of the physical fields
    end
    
    properties(Hidden=true)
        draw_h  % figure handles
    end
    
    properties(SetAccess=protected)
        f_extQ      % external field
        obc_file    % the open boundary file
        out_file    % the output file
    end
    properties
        mesh    % mesh object
        f_Q     % the physical field values (each pages)
    end
    %% 虚函数
    methods(Abstract)
        draw(obj, field) % draw the field
        [ f_Q ] = init(obj) % 初始化
    end
    
    methods(Abstract, Access=protected)
        [ E ] = flux_term( obj, f_Q ) % 计算体积积分通量项 F
        [ dflux ] = surf_term( obj, f_Q ) % 计算边界积分通量差值 (Fn - Fn*)
        [ rhs ] = rhs_term(obj, f_Q ) % 计算右端项
    end
    
    methods(Access=protected)
        function [ f_ext ] = ext_func(obj, time) % 解析解
        end% func
    end
    
    %% public
    methods
        function obj = phys(mesh)
            obj.mesh = mesh;
            obj.f_Q = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
            obj.f_extQ = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
        end% func
    
        function err = norm_err2(obj, time)
        % NORM_ERR2 calculate the L2 norm error from the exact solutions.
        % Warrning: The methods calls the public method - ext_func
        %   f_ext = ext_func(obj, t)
        % which returns the exact solution of each nodes at time 't'.
            f_ext = ext_func(obj, time); % get the exact solution
            err = zeros(obj.Nfield, 1);
            f_abs = obj.f_Q - f_ext;
            area = sum(obj.mesh.vol);
            for fld = 1:obj.Nfield
                temp = f_abs(:,:,fld).*f_abs(:,:,fld);
                err(fld) = sqrt( sum( ...
                    obj.mesh.cell_mean(temp).*obj.mesh.vol ) )./area;
            end
        end% func
        
        function err = norm_err1(obj, time)
            % 计算2范数误差。
            % 警告，调用此函数时需要首先定义精确解函数：
            %   f_ext = ext_func(obj, time)
            % ext_func 根据输入时间返回各节点精确解。
            f_ext = ext_func(obj, time);
            err = zeros(obj.Nfield, 1);
            f_abs = obj.f_Q - f_ext;
            area = sum(obj.mesh.vol);
            for fld = 1:obj.Nfield
                temp = abs( f_abs(:,:,fld) );
                err(fld) = sum( ...
                    obj.mesh.cell_mean(temp).*obj.mesh.vol )./area;
            end
        end% func
        
        function err = norm_errInf(obj, time)
            % 计算最大范数误差。
            % 警告，调用此函数时需要首先定义精确解函数：
            %   f_ext = ext_func(obj, time)
            % ext_func 根据输入时间返回各节点精确解。
            f_ext = ext_func(obj, time);
            err = zeros(obj.Nfield, 1);
            f_abs = obj.f_Q - f_ext;
            for fld = 1:obj.Nfield
                temp = abs( f_abs(:,:,fld) );
                err(fld) = max(max( temp ));
            end
        end% func
    end% methods
    
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
            vert_extQ = obj.obc_file.get_extQ(stime);
            vertlist = obj.obc_file.vert;
            for fld = 1:obj.Nfield % map vertex values to nodes
                vert_Q = zeros(obj.mesh.Nv, 1);
                vert_Q( vertlist ) = vert_extQ(:, fld);
                obj.f_extQ(:,:,fld) = obj.mesh.proj_vert2node(vert_Q);
            end
            
        end% func
    end
    
end