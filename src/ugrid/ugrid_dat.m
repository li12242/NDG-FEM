classdef ugrid_dat < handle
    %UGRID_SET This class represent a dat set in ugrid

    properties
        set_source   % the associated set
        dim   % num of variables per element
        type  % datatype
        data  % input data values
        name  % name of set variable
    end
    
    methods
        function obj = ugrid_set(set_source, dim, type, data, name)
            %UGRID Construct an instance of this class
            obj.set_source = set_source;
            obj.dim = dim;
            obj.type = type;
            obj.data = data;
            obj.name = name;
        end
    end

end % class