classdef obc_file < ndg_utility.nc.nc_file
    %OBC_FILE Summary of this class goes here
    %   Detailed explanation goes here
    
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
        
        function vertExt = interp(obj, Nfield, step, coef)
            tmp1 = netcdf.getVar(obj.ncid, obj.f_extID, ...
                [0, 0, step(1)-1], [obj.Nv, Nfield, 1]);
            tmp2 = netcdf.getVar(obj.ncid, obj.f_extID, ...
                [0, 0, step(2)-1], [obj.Nv, Nfield, 1]);
            vertExt = tmp1.*coef(1) + tmp2.*coef(2); 
        end
    end
    
    %% public methods
    methods
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
        
        function vertQ = get_extQ(obj, Nfield, stime)
            % Get the external value
            [step, coef] = coef_parameter_update(obj, stime);
            vertQ = interp(obj, Nfield, step, coef);
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
