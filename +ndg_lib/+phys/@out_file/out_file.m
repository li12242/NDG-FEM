classdef out_file < ndg_utility.nc.nc_file
    %OUT_FILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt % output time interval
        counter = 0; % counter of out step
    end
    
    methods
        declare_file(obj) % create netcdf file with contained variables
    end
    
    methods
        function obj = out_file(name, Np, K, Nfield, dt)
            % Construct empty nc_file object
            obj.name = name;
            obj.dt = dt;
            
            % define dimensions
            np = ndg_utility.nc.nc_dim('Np', Np);
            ne = ndg_utility.nc.nc_dim('K', K);
            nf = ndg_utility.nc.nc_dim('Nfield', Nfield);
            nt = ndg_utility.nc.nc_dim('time', 0);
            obj.dims = [np, ne, nf, nt];
            % define variables
            fq = ndg_utility.nc.nc_var('f_Q', obj.dims, 'float');
            ti = ndg_utility.nc.nc_var('time', nt, 'float');
            obj.vars = [fq, ti];
        end
        
        function delete(obj) % 析构函数
            if(obj.isopen) % if netcdf file is still open
                netcdf.close(obj.ncid);
                obj.isopen = false;
                obj.counter = 0;
            end
        end% func
        
        function close(obj) % close output file
            if(obj.isopen) % if netcdf file is still open
                netcdf.close(obj.ncid);
                obj.isopen = false;
                obj.counter = 0;
            end
        end
        
        function dim = add_dims(obj, name, len)
            % add dimensions
            dim = ndg_utility.nc.nc_dim(name, len);
            obj.dims = [obj.dims, dim];
        end
        
        function obj = add_vars(obj, name, dims, type)
            % add variables
            var = ndg_utility.nc.nc_var(name, dims, type);
            obj.vars = [obj.vars, var];
        end
        
        function put_var(obj, time, f_Q)
            % 将变量 f_Q 与 time 写入结果文件内
            if (time > obj.counter*obj.dt)
                netcdf.putVar(obj.ncid, obj.vars(1).id, [0, 0, obj.counter],...
                    [obj.vars(1).dims(1:2).len, 1], f_Q);
                netcdf.putVar(obj.ncid, obj.vars(2).id, obj.counter, 1, time);
                obj.counter = obj.counter + 1;
            end% if
        end
    end% methods
    
end

