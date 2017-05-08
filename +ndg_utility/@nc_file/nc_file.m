classdef nc_file < matlab.mixin.SetGet
    %NC_FILE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = private)
        isopen = false  % indicator for wheather the file is open
        ncid    % ID of ncfile
        ncdims  % array of dimensions in NetCDF file
        ncvars  % array of variables in NetCDF file
        ncname  % file name of NetCDF file
    end
    
    %% protected methods
    methods(Access=protected)
        function [] = read_file(obj, varargin)
            % Get dim and var from netcdf file
            % Usages: 
            %   file.read_file(''); % read properties from input file 
            %   file.read_file(); % equivalent to file.real_file(file.ncname)
            %   
        end
    end
    
    %% public methods
    methods
        function obj = nc_file()
            % Constructor
        end% func
        
        function obj = set.ncdims(obj, dimarray)
        end
        
        function obj = set.ncvars(obj, vararray)
        end
        
        function declare_file(obj, name)
            % Create a new NetCDF file and define with contained variables
        end
        
        function obj = set_file(obj, filename)
            % Set the object attach with a new NetCDF file
            if( exist(filename, 'file') ~= 2 ) % check file exist
                error(['Could not find file: ', filename]);
            end
            obj.ncname = filename;
            if( obj.isopen ) % if the old file is still open
                netcdf.close(obj.ncid);
            end
            obj.ncid = netcdf.open(filename); % open new file
            obj.isopen = true;
            
        end
    end
    
end

