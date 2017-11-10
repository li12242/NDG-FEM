classdef conv1d < ndg_lib.phys.phys1d
    %CONV1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        Nfield = 1
    end
    
    properties(SetAccess = protected)
        cfl     % CFL 数
        ftime   % 计算终止时间
        dt      % 计算时间步长
        u       % 速度场
    end
    
    %% 虚函数
    methods(Abstract)
        [ spe ] = character_len(obj, f_Q)
    end
    
    %% function - matlab
    methods(Access=protected) % private 
        [ E ] = flux_term( obj, f_Q ) % get the flux terms
        [ dflux ] = surf_term( obj, f_Q ) % get flux deviation
        [ rhs ] = rhs_term(obj, f_Q ) % get the r.h.s term
    end
    
    %% function - mat-mex
    methods(Access=protected) % private
    end
    methods % mex solver
    end
    
    %% private function - mat-gpu
    methods(Access=protected) % private
    end
    methods % gpu solver
    end
    
    %% private function - cuda-gpu
    methods(Access=protected) % private
    end
    methods % cuda solver
    end
    
    %% 公共函数
    methods  
        
        function obj = conv1d(mesh)
            obj = obj@ndg_lib.phys.phys1d(mesh);
        end
        
        function create_obc_file(obj, filename, dt, ftime)
            time = (0:dt:ftime);
            [minvx, minId] = min(obj.mesh.vx); % 左侧顶点
            [maxvx, maxId] = max(obj.mesh.vx); % 右侧顶点
            vx = [minvx, maxvx];
            fun = @(t) sin( (vx/obj.u(1) - t)/ftime*2*pi ); % 外部值函数
            Nt = numel(time);
            % dimensions
            flddim = ndg_utility.nc.nc_dim('Nfield', obj.Nfield);
            vdim   = ndg_utility.nc.nc_dim('Nv', 2);
            tdim   = ndg_utility.nc.nc_dim('time', Nt);
            % variables
            vertID = ndg_utility.nc.nc_var('vert', vdim, 'int');
            timeID = ndg_utility.nc.nc_var('time', tdim, 'float');
            variID = ndg_utility.nc.nc_var('f_extQ', [vdim, flddim, tdim], 'double');
            % file
            ncid = netcdf.create(filename,'CLOBBER');
            % define dims
            flddim.define_in_ncfile(ncid);
            vdim.define_in_ncfile(ncid);
            tdim.define_in_ncfile(ncid);
            % define variables
            vertID.define_in_ncfile(ncid);
            timeID.define_in_ncfile(ncid);
            variID.define_in_ncfile(ncid);
            
            netcdf.endDef(ncid);
            % put variables
            netcdf.putVar(ncid, vertID.id, [minId, maxId]);
            netcdf.putVar(ncid, timeID.id, time);
            for t = 1:Nt
                netcdf.putVar(ncid, variID.id, [0, 0, t-1], ...
                    [2, obj.Nfield, 1], fun(time(t)));
            end
            netcdf.close(ncid);
        end
    end% methods
    
end

