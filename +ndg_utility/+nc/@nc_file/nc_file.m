classdef nc_file < matlab.mixin.SetGet
    %NC_FILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = protected)
        isopen = false  % indicator for wheather the file is open
        ncid    % ID of ncfile
        dims  % array of dimensions in NetCDF file
        vars  % array of variables in NetCDF file
        name  % file name of NetCDF file
    end
    
    %% protected methods
    methods(Access=protected)
        function [obj] = read_file(obj)
            % Get dim and var from a new netcdf file
            % Usages: 
            %   file.read_file(''); % read properties from input file 
            obj.ncid = netcdf.open(obj.name, 'NOWRITE');
            
            dimids = netcdf.inqDimIDs(obj.ncid);
            dim = repmat(ndg_utility.nc.nc_dim(), numel(dimids), 1);
            for n = 1:numel(dimids)
                dim(n) = ndg_utility.nc.nc_dim();
                read_from_file(dim(n), obj.ncid, dimids(n));
            end
            obj.dims = dim;
            
            varids = netcdf.inqVarIDs(obj.ncid);
            var = repmat(ndg_utility.nc.nc_var(), numel(varids), 1);
            for n = 1:numel(varids)
                var(n) = ndg_utility.nc.nc_var();
                read_from_file(var(n), obj.ncid, varids(n), dim);
            end
            obj.vars = var;
            obj.isopen = true;
        end
    end
    
    %% public methods
    methods
        
        function obj = nc_file(varargin)
            % create nc_file object
            % Usages:
            %   file = nc_file(); % create an empty nc_file object 
            %   file = nc_file('test.nc', [x, time], [var]); % create a
            %   nc_file object with contains.
            switch nargin
                case 0
                    return;
                case 3
                    obj.name = varargin{1};
                    obj.dims = varargin{2};
                    obj.vars = varargin{3};
                otherwise
                    error(['Input variable number', num2str(nargin), 'is incorrect!']);
            end
        end% func
        
        function delete(obj)
            if(obj.isopen) % if netcdf file is still open
                netcdf.close(obj.ncid);
                obj.isopen = false;
            end
        end% func
        
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
        end
    end
    
end

