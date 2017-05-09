classdef nc_var < matlab.mixin.SetGet
    %NC_VAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        dims        % nc_dim array
        name        % variable name
        id          % variable ID
        type        % variable types
    end% properties
    
    %% private methods
    methods(Access=private)
        function set_type(obj, type)
            % set variable types as double, float, short or int.
            switch type
                case 'double'
                    obj.type = type;
                case 'float'
                    obj.type = type;
                case 'short'
                    obj.type = type;
                case 'int'
                    obj.type = type;
                otherwise
                    error(['Unknown variable type: ', type, '\n']);
            end% switch
        end% func
        
        function set_dims(obj, dims)
            % set dimension arrays
            ndim = numel(dims);
            if ndim > 1  % has more than one dims, only the length 
                % of the last dim could be unlimited
                if any( dims( 1:(ndim-1) ).len < 0 )
                    error(['The length of dimension ', ...
                        'should be positive intger!'])
                elseif any( dims( 1:(ndim-1) ).len == 0 )
                    error(['The length of dimension ', ...
                        'should not be unlimited!'])
                end% if
            else
                % only has one dimension
                if dims.len < 0
                    error(['The length of dimension ', dims.name, ...
                        'should be positive intger!'])
                end% if
            end% if
            obj.dims = dims;
        end% func
    end% methods
    
    %% public methods
    methods
        function obj = nc_var(varargin)
            switch nargin
                case 1
                    obj.name = varargin{1};
                case 3
                    obj.name = varargin{1};
                    obj.set_dims(varargin{2});
                    obj.set_type(varargin{3});
            end
        end
        
        function read_from_file(obj, ncid, varid)
            obj.id = varid;
            [obj.name,obj.type,obj.dims,~] = netcdf.inqVar(ncid, varid);
        end
        
        function define_in_ncfile(obj, ncid)
            % define variable in NetCDF file
            obj.id = netcdf.defVar(ncid, obj.name, obj.type, ...
                [obj.dims(1:end).id]);
        end% func
    end
    
end

