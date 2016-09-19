classdef StdCell < handle
    properties
        dim     % element dimension
        np      % number of points in each element
        bclist  % boundary node index
        verlist % vertex node list
        type    % name of element type
    end% properties
    
    methods
        function obj = StdCell(eleType, order)
            obj.type = eleType;
            switch eleType
                case 'line'
                    obj.dim = 1;
                    shape   = StdRegions.Line(order);
                    obj.np  = shape.nNode;
                    obj.bclist  = shape.getFaceListToNodeList;
                    obj.verlist = shape.getVertexNodeList; 
                case 'tri'
                    obj.dim = 2;
                    shape   = StdRegions.Triangle(order);
                    obj.np  = shape.nNode;
                    obj.bclist = shape.getFaceListToNodeList;
                    obj.verlist = shape.getVertexNodeList; 
                case 'quad'
                    obj.dim = 2;
                    shape   = StdRegions.Quad(order);
                    obj.np  = shape.nNode;
                    obj.bclist = shape.getFaceListToNodeList;
                    obj.verlist = shape.getVertexNodeList; 
                otherwise
                    error(['The input element type must be one of "line",',...
                        ' "tri" or "quad"'])
            end% swithch
        end% func
    end% methods
    
end% classdef