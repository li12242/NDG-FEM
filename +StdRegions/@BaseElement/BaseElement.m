classdef BaseElement < handle
    properties(SetAccess=protected)
        nDim        % the dimension of element
        nVertice    % the number of vertices
        nOrder      % the order of degrees
        nFace       % the number of edges
        nNode       % the number of nodes
        sName       % the name of element
    end
    
    methods
        function obj = BaseElement(dim, vertice, order, nface)
        % nNode sName is defined by subclass
            obj.nDim = dim;
            obj.nVertice = vertice;
            obj.nOrder = order;
            obj.nFace = nface;
        end
    end
end