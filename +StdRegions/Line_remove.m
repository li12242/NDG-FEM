classdef Line < StdRegions.LineBasic
    % Standard Element - Line
    % Line properties(Inherit): 
    %   nDim        - the dimension of element 
    %   nFace       - the number of edges 
    %   nNode       - the number of nodes
    %   nOrder      - the order of degrees
    %   nVertice    - the number of vertices
    %   r           - node coordinate, size [nNode x 1]
    %   Dr          - derivative matrix 
    %   M           - Mass Matrix
    %   invM        - inverse of Mass Matrix
    %   VandMatrix	- Vandmonde Matrix
    %   sName       - the name of element
    % Line properties:
    %   nFaceNode	- the number of nodes on faces
    %   Mef	 - face integral mass matrix of all nodes 
    %   Mes	 - face integral mass matrix of face nodes 
    %
    % Line methods:
    %   getNodeListAtFace        - return the number of nodes & node list
    %   getFaceListAtFace        - return the face list at spicific face
    %   getReorderFaceListAtFace - return the reorder face list of the spicific face
    %   getFaceListToNodeList    - return the node list of face node, from
    %   left to right
    %   getFaceGeometric         - return face normal vector & surface jacobi factor
    %   getEleGeometric          - return node Coordinate & dr/dx & jacobi factor
    %
%% properties public
% 
    properties
        nFaceNode               % nFace x nNode (at face element)
        Mes                     % Mass Matrix of Edge, integral only of face nodes, M_{E,small}
        Mef                     % Mass Matrix of Edge, integral of all nodes, M_{E,full}
        LIFT                    % lift matrix, inv(M)*Mes
        Fmask
    end% public properties
%% properties private

    properties(SetAccess=private, GetAccess=private)
        nPerBoundaryNode   % number of nodes at every boundary, size [nFace x 1] 
    end
    methods
%% function line
% 
        function obj = Line(nOrder)
            % construction function
            obj = obj@StdRegions.LineBasic(nOrder);
            PointFaceShape = StdRegions.Point();
            % 2 faces of point type
            obj.nPerBoundaryNode = [PointFaceShape.nNode; PointFaceShape.nNode];
            % nFaceNode = nFace x nNode
            obj.nFaceNode = sum(obj.nPerBoundaryNode);
            obj.Mes = getSmallFaceMassMatrix(obj, PointFaceShape);
            obj.Mef = getFullFaceMassMatrix(obj);
            obj.LIFT = (obj.VandMatrix*(obj.VandMatrix)')*obj.Mes;
            
            % Fmask
            obj.Fmask = zeros(obj.nFace, PointFaceShape.nNode);
            for i = 1:obj.nFace
                [~, t] = obj.getNodeListAtFace(i);
                obj.Fmask(i, :) = t';
            end% for
        end% func
%% function getNodeListAtFace
% get the spicific face local face list
% return the number of nodes & node indicator, at spicific face
% 
        function  [n, nodelist] = getNodeListAtFace(obj, iface)
            
            facelist = obj.getFaceListAtFace(iface);
            faceListToNodelist = getFaceListToNodeList(obj);
            nodelist = faceListToNodelist(facelist);
            n = obj.nPerBoundaryNode(iface);
        end% func
%% function getFaceListAtFace
% 
% 
        function facelist = getFaceListAtFace(obj, iface) 
            % return the face list at spicific face
            if (iface ==1)
                contour = 0; % the sum of boundary node num before iface
            else 
                contour = sum( obj.nPerBoundaryNode(1:iface-1) );
            end
            % the number of nodes at the ibth boundary is 'nPerBoundaryNode(ib)'
            facelist = contour+1:contour+obj.nPerBoundaryNode(iface);
        end %function getLocalFaceList
%% function getReorderFaceListAtFace
% 
% 
        function facelist = getReorderFaceListAtFace(obj, iface, ~)
            % return the reorder face list of the spicific face
            % with the vertice order reformed
            facelist = getFaceListAtFace(obj, iface);
        end% function
        
%%  function getFaceListToNodeList
% 
% 
        function nodelist = getFaceListToNodeList(obj)
            % return the node list of face node
            % size [nFaceNode x 1]
            nodelist = zeros(obj.nFaceNode, 1);
            nodelist(1) = 1;
            nodelist(2) = obj.nOrder+1;
            nodelist = int16(nodelist);
        end %function
%% function getFaceGeometric
% 
% 
        function [nx, sJ] = getFaceGeometric(~, x)
            % get Face Normal vector & surface jacobi factor
            % Input:    x  - node coordinate, size [nNode, nElement]
            % Output:   nx - outward normal vector
            %           sJ - surface edge jacobi coefficient
            nx = zeros(2, size(x,2));
            % Define outward normals
            nx(1, :) = -1.0; nx(2, :) = 1.0; sJ = ones(size(nx));
        end% function
%% function getEleGeometric  
% 
% 
        function [x, rx, J] = getEleGeometric(obj, vx)
            % get element Geometric Factor
            % Input:    vx  - Vertic Coordinate, size [2(nVertice) x nElement]
            % Output:   x   - node coordinate
            %           rx  - dr/dx at nodes
            %           J   - jacobi factor
            assert(size(vx,1)==2, 'transferToPhysic: input vx faults')
            x = 0.5*((1-obj.r)*vx(1,:) + (obj.r+1)*vx(2,:));
            xr  = obj.Dr*x; J = xr; rx = 1./J;
        end% function
    end %methods public
    methods(Hidden)
%% function getSmallFaceMassMatrix
% 
% 
        function FaceMassMatrixSmall = getSmallFaceMassMatrix(obj, PointShape)
            % getSmallFaceMassMatrix
            % $M^f_{i,j} = \int_{\partial \Omega}l_i l_{fj} ds$
            % size [nNode x nFaceNode]

            % allocate the Face Mass Matrix 
            FaceMassMatrixSmall = zeros(obj.nNode,obj.nFaceNode);
            for ib = 1:obj.nFace
                % iFaceList: the No. of ibth boundary local face node list. 
                iFaceList = obj.getFaceListAtFace(ib);
                [~,iFaceNodeList] = obj.getNodeListAtFace(ib);
                FaceMassMatrixSmall(iFaceNodeList,iFaceList) = PointShape.M;
            end
        end %function getSmallFaceMassMatrix
%% function getFullFaceMassMatrix
%  
% 
        function FaceMassMatrixFull = getFullFaceMassMatrix(obj)
            % getFullFaceMassMatrix
            % $M^f_{i,j} = \int_{\partial \Omega}l_i l_j ds$, size [nNode x nNode]
            
            % allocate the Face Mass Matrix 
            FaceMassMatrixFull = zeros(obj.nNode, obj.nNode);
            for ib = 1:obj.nFace
                [~,nodelist] = obj.getNodeListAtFace(ib);
                facelist = obj.getFaceListAtFace(ib);
                FaceMassMatrixFull(nodelist, nodelist) ...
                    = FaceMassMatrixFull(nodelist, nodelist) + ...
                    obj.Mes(nodelist,facelist);
            end
        end %function getFullFaceMassMatrix
    end% methods private
end %classdef Line


