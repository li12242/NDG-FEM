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
%% properties private
    properties(SetAccess=private, GetAccess=private)
        nPerBoundaryNode   % [nFace x 1] number of nodes at every boundary
    end% properties private
    
    methods
        function obj = Quad(nOrder)
            obj = obj@StdRegions.QuadBasic(nOrder);
            
            % face element
            LineFaceShape = StdRegions.LineBasic(nOrder);
            % 4 faces of line type
            obj.nPerBoundaryNode = [LineFaceShape.nNode; LineFaceShape.nNode; ...
                LineFaceShape.nNode; LineFaceShape.nNode;];
            % nFaceNode = nFace x nNode
            obj.nFaceNode = sum(obj.nPerBoundaryNode);
            obj.Mes = getSmallFaceMassMatrix(obj, LineFaceShape);
            obj.Mef = getFullFaceMassMatrix(obj);
        end% func
        
        function [x, y, rx, sx, ry, sy, J] = getEleGeometric(obj, VX, VY)
            % get element Geometric Factor
            % Input:    vx - Vertic Coordinate, size [3(nVertice) x nElement]
            %           vy - Vertic Coordinate, size [3(nVertice) x nElement]
            % Output:   x  - node coordinate
            %           rx - dr/dx at nodes
            %           J  - jacobi factor
            assert((size(VX,1)==4 && size(VY,1)==4), 'transferToPhysic: input vertice faults')
            
            % coordinate mapping
            x = (1-obj.r).*(1-obj.s)./4*VX(1,:) + (1+obj.r).*(1-obj.s)./4*VX(2,:)...
                +(1+obj.r).*(1+obj.s)./4*VX(3,:) +(1-obj.r).*(1+obj.s)./4*VX(4,:);
            y = (1-obj.r).*(1-obj.s)./4*VY(1,:) + (1+obj.r).*(1-obj.s)./4*VY(2,:)...
                +(1+obj.r).*(1+obj.s)./4*VY(3,:) +(1-obj.r).*(1+obj.s)./4*VY(4,:);
            
            xr = obj.Dr*x; xs = obj.Ds*x; yr = obj.Dr*y; ys = obj.Ds*y; 
            J = -xs.*yr + xr.*ys;
            rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
        end
        
        function [nx, ny, sJ] = getFaceGeometric(obj, x, y)
            % get Face Normal vector & surface jacobi factor
            % Input:    x - node coordinate, size [nNode, nElement]
            %           y - node coordinate, size [nNode, nElement]
            % Output:   nx - outward vector
            %           ny - outward vector
            %           sJ - face jacobi factor
            K = size(x, 2);
            xr = obj.Dr*x; yr = obj.Dr*y; xs = obj.Ds*x; ys = obj.Ds*y;
            % interpolate geometric factors to face nodes
            Fmask = obj.getFaceListToNodeList;
            fxr = xr(Fmask, :); fxs = xs(Fmask, :); 
            fyr = yr(Fmask, :); fys = ys(Fmask, :);
            % build normals
            Nfp = obj.nFaceNode/4;
            Nf = 4;
            nx = zeros(Nf*Nfp, K); ny = zeros(Nf*Nfp, K);
            fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; fid3 = (2*Nfp+1:3*Nfp)';
            fid4 = (3*Nfp+1:4*Nfp)';
            % face 1
            nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);
            % face 2
            nx(fid2, :) =  fys(fid2, :); ny(fid2, :) = -fxs(fid2, :);
            % face 3
            nx(fid3, :) = -fyr(fid3, :); ny(fid3, :) =  fxr(fid3, :);
            % face 4
            nx(fid4, :) = -fys(fid4, :); ny(fid4, :) =  fxs(fid4, :);
            % normalise
            sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
        end
        
        function facelist = getFaceListAtFace(obj,iface)
        %return the face list at spicific face
        n = obj.nOrder;
        np = n+1;
        facelist = (iface - 1)*np+1:iface*np;
        end% func
        
        function [n, nodelist] = getNodeListAtFace(obj,iface)
            facelist = obj.getFaceListAtFace(iface);
            NodeListAtFace = getFaceListToNodeList(obj);
            nodelist = NodeListAtFace(facelist);
            n = obj.nPerBoundaryNode(iface);
        end% function
        
        function nodelist = getFaceListToNodeList(obj)
            nf = obj.nFace; % No of faces
            np = obj.nOrder+1; % No of points at one face
            nodelist = zeros(nf*np, 1);
            
            % node arrangement matrix
            % 
            indMatrix = flip(reshape(1:np^2, np, np)');
            
            % face I, bottom
            nodelist(1:np) = indMatrix(np, :);
            % face II, right
            nodelist(np+1:2*np) = flip(indMatrix(:, np));
            % face III, up
            nodelist(2*np+1:3*np) = flip(indMatrix(1, :));
            % face IV, left
            nodelist(3*np+1:4*np) = indMatrix(:, 1);
        end% func
    end
    
    methods(Hidden)
        function FaceMassMatrixSmall = getSmallFaceMassMatrix(obj, LineShape)
            % getSmallFaceMassMatrix
            % $M^f_{i,j} = \int_{\partial \Omega}l_i l_{fj} ds$, 
            % size [nNode x nFaceNode]

            % allocate the Face Mass Matrix 
            FaceMassMatrixSmall = zeros(obj.nNode, obj.nFaceNode);
            for ib = 1:obj.nFace
                % iFaceList: the No. of ibth boundary local face node list. 
                iFaceList = obj.getFaceListAtFace(ib);
                [~,iFaceNodeList] = obj.getNodeListAtFace(ib);
                FaceMassMatrixSmall(iFaceNodeList,iFaceList) = LineShape.M;
            end
        end %function getSmallFaceMassMatrix
        
        
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
    end% methods
end% class