classdef NcDim < handle
    properties
        name    % dimension name
        len     % legnth
        id      % id in nctcdf file
    end% properties
    
    methods
        
        %% dimobj
        % Define dimensions for NetCDF file
        % Input:
        %   name - name of dimension
        %   len  - dimension length (0 for unlimited dimensions)
        %
        % Usages
        %   time = Utilities.netcdf.dimobj('time', 0)
        function obj = NcDim(name, len)
            obj.name = name;
            if len >= 0
                obj.len = len;
            else
                error(['The length of dimension ', name, 'should be positive intger!\n'])
            end
        end% func
        
        function obj = setID(obj, id)
            % Set the dimension ID
            obj.id = id;
        end% func
        
    end% methods
end% class