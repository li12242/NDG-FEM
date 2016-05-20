classdef varobj < handle
    properties
        dimArray
        name
        id % ID in the netcdf file
        type % variable types
    end% properties
    
    methods
        
        function obj = varobj(name, dimArray, type)
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
            obj.type = type;
        end% func
        
        function obj = setID(obj, id)
            obj.id = id;
        end% func
    end% methods
end% class