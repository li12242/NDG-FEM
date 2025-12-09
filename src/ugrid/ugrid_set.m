classdef ugrid_set < handle
    %UGRID_SET This class represent a set in ugrid
    
    properties
        size  % length of set
        name  % name of set variable
    end
    
    methods
        function obj = ugrid_set(size, name)
            %UGRID Construct an instance of this class
            obj.size = size;
            obj.name = name;
        end
    end
end

