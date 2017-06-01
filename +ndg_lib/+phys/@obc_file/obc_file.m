classdef obc_file < ndg_utility.nc.nc_file
    %OBC_FILE 开边界文件类
    %   开边界文件默认为 NetCDF 格式，通过时间线性插值读取任意时刻开边界数据。
    %   对象包含的方法有:
    %       obc_file - 构造函数，包括无输入参数和输出开边界文件名两种调用方式;
    %       get_extQ - 获取开边界文件内外部值数据，其中外部值为顶点处二维数据，
    %           其各个维度分别为 [vert, field];
    %       set_file - 指定开边界文件，此方法用于辅助生成 obc_file 对象为空的情况;
    %       make_obc_file - 通过输入顶点，时间，以及外部值的方法生成开边界 NetCDF
    %           文件，其中外部值为三维数据，其各个维度分别为 [vert, field, time];
    
    properties(Constant)
        Nstep = 2 % use adjacent two steps to interpolate
    end
    
    properties(SetAccess=protected)
        Ntime   % 开边界文件时间步
        time    % 开边界文件时间序列
        Nv      % 开边界文件顶点数
        vert    % 开边界文件顶点号,从1开始
        Nfield  % 开边界文件变量数
        f_ext
        f_extID % 变量ID
    end
    
    %% private methods
    methods(Access=private)
        function [step, coef] = coef_parameter_update(obj, stime)
            % 根据给定时间，寻找最接近的两个时间步，并给出各自插值系数.
            t = find( obj.time > stime, 1 );
            if (isempty(t)) 
                step = [obj.Ntime, obj.Ntime]; coef = [1, 0]; return; 
            end
            if (t>1)
                step = [t-1, t];
            elseif (t == 1)
                step = [t, t]; coef = [1, 0]; return;
            end
            % find coeff
            coef = abs( obj.time(step) - stime )/diff( obj.time(step) );
            coef = flip(coef);
        end% func
        
        function vertExt = interp(obj, step, coef)
            tmp1 = netcdf.getVar(obj.ncid, obj.f_extID, ...
                [0, 0, step(1)-1], [obj.Nv, obj.Nfield, 1]);
            tmp2 = netcdf.getVar(obj.ncid, obj.f_extID, ...
                [0, 0, step(2)-1], [obj.Nv, obj.Nfield, 1]);
            vertExt = tmp1.*coef(1) + tmp2.*coef(2); 
        end
    end
    
    %% public methods
    methods
        % 生成开边界 NetCDF 文件
        make_obc_file(obj, filename, Nfield, time, vert, f_extQ);
        
        function obj = obc_file(varargin)
            switch nargin
                case 0
                    return;
                case 1
                    obj.name = varargin{1};
                    obj = set_file(obj, obj.name);
            end
        end% func
        
        function vertQ = get_extQ(obj, stime)
            % 获取外部值数据
            [step, coef] = coef_parameter_update(obj, stime);
            vertQ = interp(obj, step, coef);
        end
                
        function obj = set_file(obj, filename)
            % Set the object associate with a new NetCDF file
            if( exist(filename, 'file') ~= 2 ) % check file exist
                error(['Could not find file: ', filename]);
            end
            obj.name = filename;
            if( obj.isopen ) % if the old file is still open
                netcdf.close(obj.ncid);
            end
            obj = read_file(obj); % open new file
            obj.isopen = true;
            % read Nfield
            nfldID = netcdf.inqDimID(obj.ncid, 'Nfield');
            [~, obj.Nfield] = netcdf.inqDim(obj.ncid, nfldID);
            % read time
            timeID = netcdf.inqVarID(obj.ncid, 'time'); % get time ID
            obj.time = netcdf.getVar(obj.ncid, timeID); % read time
            obj.Ntime = numel(obj.time);
            % read vertex
            vertID = netcdf.inqVarID(obj.ncid, 'vert'); % get vertex ID
            obj.vert = netcdf.getVar(obj.ncid, vertID); % read vertex
            obj.Nv = numel(obj.vert);
            % read field ID
            obj.f_extID = netcdf.inqVarID(obj.ncid, 'f_extQ'); % get var ID
        end
    end% methods
    
end
