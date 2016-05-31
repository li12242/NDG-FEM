classdef QuadBasic < StdRegions.BaseElement
    % Stdandard Element Quadrilatral Basic
    % QuadBasic properties(Inherit):
    %   nDim        - the dimension of element 
    %   nFace       - the number of edges 
    %   nNode       - the number of nodes 
    %   nOrder      - the order of degrees 
    %   nVertice	- the number of vertices 
    %   sName       - the name of element
    % QuadBasic properties:
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
        
        Dr  % Derivative Matrix of r
        Ds  % Derivative Matrix of s
        Drw % Derivative Matrix of r in weak form 
        Dsw % Derivative Matrix of s in weak form 
        r   % point coordinate of dim 1, [nNode x 1] 
        s   % point coordinate of dim 2, [nNode x 1] 
        VandMatrix  % Vandmonde Matrix
    end% properties
    
    methods
        function obj = QuadBasic(nOrder)
            import StdRegions.Quad.*
            
            dim = 2; vertice = 4; face = 4;
            obj = obj@StdRegions.BaseElement(dim,vertice,nOrder,face);
            obj.sName = 'Quad';
            obj.nNode = (nOrder+1)^2;
            
            % get coordinate of LGL points on reference element
            [obj.r, obj.s] = getCoor(nOrder);
            
            % get Vandermonde matrix
            obj.VandMatrix = getVandMatrix(obj.nOrder, obj.r, obj.s);
            invV = inv(obj.VandMatrix);
            obj.M = invV'*invV;
            
            
            [obj.Dr, obj.Ds, obj.Drw, obj.Dsw] = ...
                getDeriMatrix(obj.nOrder, obj.r, obj.s ,obj.VandMatrix);
        end% func
    end
end% classdef