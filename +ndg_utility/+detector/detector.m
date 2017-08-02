classdef detector < handle
    %DETECTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Hidden, SetAccess=protected)
        contour = 0; % 当前记录时间步;
        intMat; % 插值矩阵
    end
    
    properties(SetAccess=protected)
        mesh % 网格对象;
        Nd % 检测点个数（根据输入检测点个数确定）;
        Nt % 时间序列个数（根据时间步长与结束时间确定）;
        Nfield % 物理量个数;
        xd, yd, zd % 检测点坐标;
        kd % 检测点所在单元序号;
        rd, sd, td % 检测点所在单元局部坐标;
        dt % 时间步长;
        ftime % 结束时间
        time % 记录时间序列;
        dQ % 检测点物理量值;
    end
    
    methods(Abstract, Hidden)
        [kd, rd, sd, td] = findlocate(obj) % 寻找节点所在单元编号与单元内局部坐标
    end
    
    methods(Hidden, Access=protected)
        function intMat = interpmatrix(obj)
            % 计算每个检测点插值矩阵
            Vg = zeros(obj.Nd, obj.mesh.cell.Np);
            for n = 1:obj.mesh.cell.Np
                Vg(:, n) = obj.mesh.cell.orthogonal_func(obj.mesh.cell.N, ...
                    n, obj.rd, obj.sd, obj.td);
            end% for
            intMat = Vg/obj.mesh.cell.V;
        end% func
    end
    
    methods
        function obj = detector(mesh, xd, yd, zd, dt, ftime, Nfield)
            % 根据输入节点坐标及时间步构造检测器
            obj.mesh = mesh;
            obj.xd = xd; 
            obj.yd = yd; 
            obj.zd = zd;
            obj.Nd = numel(xd);
            obj.dt = dt; % 检测数据时间最小步长
            obj.ftime = ftime;
            obj.Nt = ceil(ftime/dt);
            obj.time = zeros(obj.Nt, 1);
            obj.Nfield = Nfield;
            
            obj.init();
            [obj.kd, obj.rd, obj.sd, obj.td] = obj.findlocate();
            obj.intMat = obj.interpmatrix();
        end% func
        
        function obj = init(obj)
            % 将检测结果初始化
            obj.contour = 0;
            obj.dQ = zeros(obj.Nd, obj.Nt, obj.Nfield);
        end% func
        
        function p_h = draw(obj, pointID, fldID)
            % 绘制检测结果
            p_h = plot(obj.time(1:obj.contour), ...
                obj.dQ(pointID, 1:obj.contour, fldID), 'b.-');
        end
        
        function collect(obj, f_Q, time)
            % 根据计算结果计算检测点插值数据
            if obj.contour > 0
                if ( (time >= obj.time(obj.contour) + obj.dt) || ...
                        (abs(time - obj.ftime) < 1e-10) )
                    obj.contour = obj.contour + 1;
                else
                    return;
                end
            else
                obj.contour = 1;
            end
            obj.time(obj.contour) = time;
            for fld = 1:obj.Nfield
                for n = 1:obj.Nd % 对每个检测点所在单元结果进行插值，获得检测点结果
                    obj.dQ(n, obj.contour, fld) = ...
                        obj.intMat(n,:)*f_Q(:,obj.kd(n),fld);
                end
            end
        end% func
    end
    
end

