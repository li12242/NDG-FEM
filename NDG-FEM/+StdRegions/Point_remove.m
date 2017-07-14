classdef Point < StdRegions.BaseElement
    properties
        M
    end% properties
    methods
        %% Point
        function obj = Point()
            dim = 0; vertice = 1; face = 0; nOrder = 0;
            obj = obj@StdRegions.BaseElement(dim,vertice,nOrder,face);
            obj.sName = 'Point';
            obj.nNode = 1;
            obj.M = 1;
        end% function
    end% methods
end% classdef