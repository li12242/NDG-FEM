classdef StdCell < handle
    properties
        dim     % element dimension
        np      % number of points in each element
        bclist  % boundary node index
        type    % name of element type
    end% properties
    
    methods
        function obj = StdCell(eleType, order)
            obj.type = eleType;
            switch eleType
                case 'line'
                    obj.dim = 1;
                    line   = StdRegions.Line(order);
                    obj.np = line.nNode;
                    obj.bclist = line.getFaceListToNodeList;
                case 'tri'
                    obj.dim = 2;
                    tri    = StdRegions.Triangle(order);
                    obj.np = tri.nNode;
                    obj.bclist = tri.getFaceListToNodeList;
                case 'quad'
                    obj.dim = 2;
                    quad   = StdRegions.Quad(order);
                    obj.np = quad.nNode;
                    obj.bclist = quad.getFaceListToNodeList;
                otherwise
                    error(['The input element type must be one of "line",',...
                        ' "tri" or "quad"'])
            end% swithch
        end% func
    end% methods
    
end% classdef