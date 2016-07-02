classdef NcVar < handle
    properties
        dimArray    % variable dimension array
        name        % variable name
        id          % variable ID
        type        % variable types
    end% properties
    
    methods
        %% NcVar
        % Define variables for NetCDF output.
        %
        % Usages:
        %
        %   time = Utilities.Netcdf.NcDim('time', 0);
        %   node = Utilities.Netcdf.NcDim('node', 10);
        %   x    = Utilities.Netcdf.NcVar('x'   , [node time], 'double');
        %
        function obj = NcVar(name, dimArray, type)
            % variable name
            obj.name = name;
            
            % variable dimensions
            obj.dimArray = dimArray;
            
            % check dimensional
            ndim = numel(dimArray);
            if ndim > 1 
                % has more than one dims, only the length of the last dim
                % could be unlimited
                for i = 1:ndim-1
                    if obj.dimArray(i).len < 0
                        error(['The length of dimension ', obj.dimArray(i).name, ...
                        'should be positive intger!\n'])
                    elseif obj.dimArray(i).len == 0
                        error(['The length of dimension ', obj.dimArray(i).name, ...
                        'should not be unlimited!\n'])
                    end% if
                end% for
            else
                % only has one dimension
                if obj.dimArray.len < 0
                    error(['The length of dimension ', obj.dimArray.name, ...
                        'should be positive intger!\n'])
                end% if
            end% if
            
            % variable types: double, float, short, int
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
        
        %% setID
        % Set variable ID
        %
        function obj = setID(obj, id)
            obj.id = id;
        end% func
    end% methods
end% class