classdef Quad < StdRegions.QuadBasic
    % Stdandard Element Quadrilatral
    % Quad properties(Inherit):
    %   nDim        - the dimension of element
    %   nFace       - the number of edges 
    %   nNode       - the number of nodes 
    %   nOrder      - the order of degrees 
    %   nVertice	- the number of vertices 
    %   sName       - the name of element
    %   M           - Mass matrix
    %   invM        - inverse of Mass Matrix
    %   Dr          - Derivative Matrix of r
    %   Ds          - Derivative Matrix of s
    %   Drw         - Derivative Matrix of r in weak form 
    %   Dsw         - Derivative Matrix of s in weak form 
    %   r           - point coordinate of dim 1, [nNode x 1] 
    %   s           - point coordinate of dim 2, [nNode x 1] 
    % Quad properties:
    %   nFaceNode           - nFace x nNode (at face element)
    %   Mes                 - face integral mass matrix of face nodes
    %   Mef                 - face integral mass matrix of all nodes
    % Quad methods:
    %   getNodeListAtFace        - return the number of nodes & node list
    %   getFaceListAtFace        - return the face list at spicific face
    %   getReorderFaceListAtFace - return the reorder face list of the spicific face
    %   getFaceListToNodeList    - return the node list of face node
    %   getFaceGeometric         - return face normal vector & surface jacobi factor
    %   getEleGeometric          - return node Coordinate & dr/dx & jacobi factor
    % Usages:
    %   obj = Quad(nOrder)
%% properties public
    properties
        nFaceNode           % nFace x nNode (at face element)
        Mes                 % face integral mass matrix of face nodes
        Mef                 % face integral mass matrix of all nodes
    end% properties
end% class