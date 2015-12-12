classdef RegionLine < MultiRegions.Region
    % physical meshes with Line element
    % RegionLine properties(Inherit):
    %   Shape    - base element type
    %   nElement - No. of element
    %   nNode    - No. of node
    %   vmapM    - face list to node list (in mesh)
    %   vmapP    - face list to neighbour node list (in mesh)
    %   sJ       - face factor: surface Line element jacobi factor
    %   J        - element factor: element Jacobi factor
    %   EToE     - element to element
    %   EToF     - element to face
    % RegionLine properties(private):
    %   SpFToV   - face to vertices
    % RegionLine properties(public):
    %   x        - node coordinate
    %   rx       - element factor: dr/dx
    %   nx       - face factor: normal vector
    %   fScale   - face factor: sJ/J(surface node)
    % RegionLine methods: none
    % 
    properties       
        nx              % face factor: normal vector
        fScale          % face factor: sJ/J(surface node)
        rx              % element factor: dr/dx
        x               % node coordinate
    end
    
    methods        
        function obj = RegionLine(line, EToV, VX)
            % Input:    VX   - vertice coordinate, size [nVertic  x1]
            %           EToV - element to vertice, size [nElement x2]
            obj = obj@MultiRegions.Region(line, EToV);
            [obj.EToE, obj.EToF] = Connect1D(obj,EToV);
            [obj.vmapM,obj.vmapP] = BuildMap(obj,EToV, obj.EToE, obj.EToF);
            % get node coordinate
            vx = zeros(2, obj.nElement);
            vx(1,:) = VX(EToV(:,1)); vx(2,:) = VX(EToV(:,2));
            % 
            [obj.x, obj.rx, obj.J] = line.getEleGeometric(vx);
            [obj.nx, obj.sJ] = line.getFaceGeometric(obj.x);
            % face scal, used for surface integral
            obj.fScale = obj.sJ./obj.J(obj.vmapM);
        end% function
    end% methods public
    
    methods(Hidden)
        function [EToE, EToF] = Connect1D(obj, EToV)
            % Build global connectivity arrays for 1D grid based on EToV 
            Nfaces = obj.Shape.nFace;
            % Find number of elements and vertices
            K = obj.nElement; TotalFaces = Nfaces*K; Nv = obj.Shape.nVertice;
            % List of local face to local vertex connections
            vn = [1;2];
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
            SpFToF = obj.SpFToV*obj.SpFToV' - speye(TotalFaces);
            % Find complete face to face connections
            [faces1, faces2] = find(SpFToF==1);
            % Convert face global number to element and face numbers
            element1 = floor( (faces1-1)/Nfaces )  + 1;
            face1    =   mod( (faces1-1), Nfaces ) + 1;
            element2 = floor( (faces2-1)/Nfaces )  + 1;
            face2    =   mod( (faces2-1), Nfaces ) + 1;
            % Rearrange into Nelements x Nfaces sized arrays
            ind = sub2ind([K, Nfaces], element1, face1);
            EToE      = (1:K)'*ones(1,Nfaces);
            EToF      = ones(K,1)*(1:Nfaces);
            EToE(ind) = element2; EToF(ind) = face2;
        end
    
        function [vmapM,vmapP] = BuildMap(obj, EToV, EToE, EToF)
            % contour the No of face node in element
            nFaceNode = obj.Shape.nFaceNode;
            vmapM = zeros(nFaceNode, obj.nElement);
            vmapP = zeros(nFaceNode, obj.nElement);
            for ie = 1:obj.nElement
                vmapM(:,ie) = obj.Shape.getFaceListToNodeList + ...
                    (ie-1)*obj.Shape.nNode;
            end
            for ie = 1:obj.nElement
                for iface = 1:obj.Shape.nFace
                    % local face list
                    localfacelist = obj.Shape.getFaceListAtFace(iface);
                    % get next face num
                    nextEle = EToE(ie, iface);
                    nextFace = EToF(ie, iface);
                    % get vertice order
                    vn = [1; 2];
                    % get vertice order
                    thisOrder = EToV(ie, vn(iface,:));
                    nextOrder = EToV(nextEle, vn(nextFace,:));
                    % nextOrder(vorder) = thisOrder
                    vorder = Utilities.PairData(thisOrder,nextOrder);
                    nextlocalfacelist = ...
                        obj.Shape.getReorderFaceListAtFace(nextFace, vorder);
                    vmapP(localfacelist,ie) = vmapM(nextlocalfacelist, nextEle);
                end
            end
        end% function
    end %methods private
end %classdef Region1D
