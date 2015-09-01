classdef RegionTri < MultiRegions.Region
    % physical meshes with Line element
    % RegionTri properties(inherit):
    %   Shape    - base element type
    %   nElement - No. of element
    %   nNode    - No. of node
    %   vmapM    - face list to node list (in mesh)
    %   vmapP    - face list to neighbour node list (in mesh)
    %   sJ       - face factor: surface element Jacobi
    %   J        - element factor: Jacobi transfer factor
    %   EToE     - element to element
    %   EToF     - element to face
    % RegionTri properties:
    %   x       - node coordinate
    %   y       - node coordinate
    %   rx      - element factor: dr/dx
    %   sx      - element factor: ds/dx
    %   ry      - element factor: dr/dy
    %   sy      - element factor: ds/dy
    %   nx      - face factor: normal vector
    %   ny      - face factor: normal vector
    %   fScale  - face factor: sJ/J(surface node)
    % RegionTri methods: none
    properties
        x       % node coordinate
        y       % node coordinate
        rx              % element factor: dr/dx
        sx              % element factor: ds/dx
        ry              % element factor: dr/dy
        sy              % element factor: ds/dy
        nx              % face factor: normal vector
        ny              % face factor: normal vector
        fScale          % face factor: sJ/J(surface node)
    end% properties
    methods
        function obj = RegionTri(shape, EToV, VX, VY)
            % Input:    EToV - element to vertice, size [nElementx3]
            %           VX - vertice coordinate, size [nVerticx1]
            %           VY - vertice coordinate, size [nVerticx1]
            obj = obj@MultiRegions.Region(shape, EToV);
            [obj.EToE, obj.EToF] = Connect2D(obj,EToV);
            [obj.vmapM,obj.vmapP] = BuildMap(obj,EToV, obj.EToE, obj.EToF);
            
            % get node coordiante
            vx = VX(EToV'); vy = VY(EToV');
            [obj.x, obj.y, obj.rx, obj.sx, obj.ry, obj.sy, obj.J]  = shape.getEleGeometric(vx, vy);
            [obj.nx, obj.ny, obj.sJ] = shape.getFaceGeometric(obj.x, obj.y);
            
            obj.fScale = obj.sJ./obj.J(obj.vmapM);
        end% function
    end% methods
    
    methods(Hidden)
        function [vmapM,vmapP] = BuildMap(obj, EToV, EToE, EToF)
            % contour the No of face node in element
            nFaceNode = obj.Shape.nFaceNode;
            vmapM = zeros(nFaceNode, obj.nElement);
            vmapP = zeros(nFaceNode, obj.nElement);
            for ie = 1:obj.nElement
                vmapM(:,ie) = obj.Shape.getFaceListToNodeList + (ie-1)*obj.Shape.nNode;
            end
            for ie = 1:obj.nElement
                for iface = 1:obj.Shape.nFace
                    % local face list
                    localfacelist = obj.Shape.getFaceListAtFace(iface);
                    % get next face num
                    nextEle = EToE(ie, iface);
                    nextFace = EToF(ie, iface);
                    % get vertice order
                    vn = [[1,2];[2,3];[3,1]];
                    thisOrder = EToV(ie, vn(iface,:));
                    nextOrder = EToV(nextEle, vn(nextFace,:));
                    % nextOrder(vorder) = thisOrder
                    vorder = Utilities.PairData(thisOrder,nextOrder);
                    nextlocalfacelist = obj.Shape.getReorderFaceListAtFace(nextFace, vorder);
                    vmapP(localfacelist,ie) = vmapM(nextlocalfacelist, nextEle);
                end
            end
        end %function
        
        % Connect 2D
        function [EToE, EToF] = Connect2D(obj,EToV)
            Nfaces = 3;
            % Find number of elements and vertices
            K = obj.nElement; Nv = max(max(EToV));
            % Create face to node connectivity matrix
            TotalFaces = Nfaces*K;
            % List of local face to local vertex connections
            vn = [[1,2];[2,3];[1,3]];
            % Build global face to node sparse array
            obj.SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
            sk = 1;
            for k=1:K
              for face=1:Nfaces
                obj.SpFToV(sk, EToV(k, vn(face,:))) = 1;
                sk = sk+1;
              end
            end
            % Build global face to global face sparse array
            SpFToF = obj.SpFToV*obj.SpFToV' - 2*speye(TotalFaces);
            % Find complete face to face connections
            [faces1, faces2] = find(SpFToF==2);
            % Convert face global number to element and face numbers
            element1 = floor( (faces1-1)/Nfaces ) + 1;
            face1    =   mod( (faces1-1),Nfaces ) + 1;
            element2 = floor( (faces2-1)/Nfaces ) + 1;
            face2    =   mod( (faces2-1),Nfaces ) + 1;
            % Rearrange into nElements x Nfaces sized arrays
            ind = sub2ind([K, Nfaces], element1, face1);
            EToE = (1:K)'*ones(1,Nfaces); EToF = ones(K,1)*(1:Nfaces);
            EToE(ind) = element2; EToF(ind) = face2;
        end% function
    end% methods
end% classedf