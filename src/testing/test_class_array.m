classdef test_class_array < handle
    properties
        obj_array
    end
    
    methods
        function obj = test_class_array(n)
            a = array(3);
            for i =1:n
                obj.obj_array{i} = a;
            end
        end
    end %methods
    
end