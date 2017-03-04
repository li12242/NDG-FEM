classdef TriangleBasic < StdRegions.BaseElement
    % Stdandard Element Triangle Basic
    % TriangleBasic properties(Inherit):
    %   nDim        - the dimension of element 
    %   nFace       - the number of edges 
    %   nNode       - the number of nodes 
    %   nOrder      - the order of degrees 
    %   nVertice	- the number of vertices 
    %   sName       - the name of element
    % TriangleBasic properties:
    %   M           - Mass matrix
    %   invM        - inverse of Mass Matrix
    %   Dr          - Derivative Matrix of r
    %   Ds          - Derivative Matrix of s
    %   Drw         - Derivative Matrix of r in weak form 
    %   Dsw         - Derivative Matrix of s in weak form 
    %   r           - point coordinate of dim 1, [nNode x 1] 
    %   s           - point coordinate of dim 2, [nNode x 1] 
    %   VandMatrix  - Vandmonde Matrix
    properties
        M   % Mass matrix
        invM        % inverse of Mass Matrix
        Dr  % Derivative Matrix of r
        Ds  % Derivative Matrix of s
        Drw % Derivative Matrix of r in weak form 
        Dsw % Derivative Matrix of s in weak form 
        r   % point coordinate of dim 1, [nNode x 1] 
        s   % point coordinate of dim 2, [nNode x 1] 
        VandMatrix  % Vandmonde Matrix
    end% properties
    
    methods
        function obj = TriangleBasic(nOrder)
            import StdRegions.Triangle.*
            
            dim = 2; vertice = 3; face = 3;
            obj = obj@StdRegions.BaseElement(dim,vertice,nOrder,face);
            obj.sName = 'Triangle';
            obj.nNode = (nOrder+1)*(nOrder+2)/2;
            
            [obj.r, obj.s] = getCoor(nOrder);
            obj.VandMatrix = getVandMatrix(obj.nOrder, obj.r, obj.s);
            invV = inv(obj.VandMatrix);
            obj.M = invV'*invV;
            obj.invM = obj.VandMatrix*(obj.VandMatrix)';
            [obj.Dr, obj.Ds, obj.Drw, obj.Dsw] = ...
                getDeriMatrix(obj.nOrder,obj.r, obj.s ,obj.VandMatrix);
        end% function
    end% methods
end% classdef