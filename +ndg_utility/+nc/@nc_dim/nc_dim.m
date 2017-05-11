classdef nc_dim < matlab.mixin.SetGet
    %NC_DIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        name    % dimension name
        len     % positive integer, 0 for unlimited
        id      % id in nctcdf file
    end% properties
    
    %% public
    methods
        function [ obj ] = nc_dim(varargin)
            % create nc_dim object with various input
            % Usages:
            %   tdim = nc_dim('time');
            %   tdim = nc_dim('time', 0);
            switch nargin
                case 1
                    return;
                case 2
                    set_name(obj, varargin{1});
                    set_len(obj, varargin{2});
                otherwise
            end
        end% func
        
        function obj = set_name(obj, name)
            obj.name = name;
        end
        
        function obj = set_len(obj, len)
            % set dimension length
            if (len >= 0)
                obj.len = len;
            else
                error(['The length of dimension ', obj.name,...
                    'should be positive intger!\n'])
            end
        end
        
        function read_from_file(obj, ncid, id)
            obj.id = id;
            [obj.name, obj.len] = netcdf.inqDim(ncid,id);
        end
        
        function define_in_ncfile(obj, ncid)
            % define the dimension in NetCDF file
            if (obj.len > 0)
                obj.id = netcdf.defDim(ncid, obj.name, obj.len);
            else
                typeid = netcdf.getConstant('NC_UNLIMITED');
                obj.id = netcdf.defDim(ncid, obj.name, typeid);
            end
        end% func
        
        
    end% methods
    
end

