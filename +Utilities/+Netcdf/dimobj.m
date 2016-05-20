classdef dimobj < handle
    properties
        name % dimension name
        len % legnth
        id % id in nctcdf file
    end% pro
    
    methods
        
        function obj = dimobj(name, len)
            obj.name = name;
            
            if len >= 0
                obj.len = len;
            else
                error(['The length of dimension ', name, 'should be positive intger!\n'])
            end
        end% func
        
        function obj = setID(obj, id)
            obj.id = id;
        end% func
    end% methods
end% class