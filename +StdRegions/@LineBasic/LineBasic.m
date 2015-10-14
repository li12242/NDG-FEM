classdef LineBasic < StdRegions.BaseElement
    % Standard Basic Element
    % LineBasic properties(inherited): 
    %   nDim        - the dimension of element 
    %   nFace       - the number of edges 
    %   nNode       - the number of nodes
    %   nOrder      - the order of degrees
    %   nVertice    - the number of vertices
    %   sName       - the name of element
    % LineBasic properties:
    %   M           - Mass Matrix
    %   invM        - inverse of Mass Matrix
    %   Dr          - derivative matrix 
    %   r           - node coordinate, size [nNode x 1]
    %   VandMatrix	- Vandmonde Matrix
    % LineBaisc methods: none
    %
    properties
        M           % Mass Matrix
        invM        % inverse of Mass Matrix
        Dr          % derivative matrix
        r           % node coordinate, size [nNode x 1]
        VandMatrix  % Vandmonde Matrix
    end% properties
    methods
        function obj = LineBasic(nOrder)
            % construction function
            dim = 1; vertice = 2; face = 2;
            obj = obj@StdRegions.BaseElement(dim,vertice,nOrder,face);
            obj.sName = 'Line';
            obj.nNode = nOrder+1;
            % get node coordiantes
            obj.r = StdRegions.Line.getCoor(nOrder);
            % get Vamdmonde Matrix
            obj.VandMatrix = StdRegions.Line.getVandMatrix(nOrder, obj.r);
            invV = inv(obj.VandMatrix);
            % get Mass Matrix
            obj.M = invV'*invV;
            obj.invM = obj.VandMatrix*(obj.VandMatrix)';
            obj.Dr = StdRegions.Line.getDeriMatrix(obj.r);
        end% function
    end% methods
end% classdef LineBasic