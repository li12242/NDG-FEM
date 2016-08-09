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
            
            % get node coordiante
            vx = VX(EToV'); vy = VY(EToV');
            [obj.x, obj.y, obj.rx, obj.sx, obj.ry, obj.sy, obj.J]  ...
                = shape.getEleGeometric(vx, vy);
            [obj.nx, obj.ny, obj.sJ] = shape.getFaceGeometric(obj.x, obj.y);
            
            [obj.EToE, obj.EToF] = Connect2D(obj,EToV);
            [obj.vmapM,obj.vmapP] = BuildMap(obj, VX, VY, EToV, obj.EToE, obj.EToF);
            
            obj.fScale = obj.sJ./obj.J(obj.vmapM);
        end% function
    end% methods
    
    methods(Hidden)
        function [vmapM,vmapP] = BuildMap(obj, VX, VY, EToV, EToE, EToF)
            % contour the No of face node in element
            NODETOL = 1e-10;
            nFaceNode = obj.Shape.nFaceNode;
            vmapM = zeros(nFaceNode, obj.nElement);
            vmapP = zeros(nFaceNode, obj.nElement);
            for k1 = 1:obj.nElement
                vmapM(:,k1) = obj.Shape.getFaceListToNodeList + (k1-1)*obj.Shape.nNode;
            end
            
            one = ones(1, obj.Shape.nFaceNode./obj.Shape.nFace);
            for k1 = 1:obj.nElement
                for f1 = 1:obj.Shape.nFace

                    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
                    
                    flist1 = obj.Shape.getFaceListAtFace(f1);
                    flist2 = obj.Shape.getFaceListAtFace(f2);

                    % reference length of edge
                    v1 = EToV(k1,f1); v2 = EToV(k1, 1+mod(f1, obj.Shape.nFace));
                    refd = sqrt( (VX(v1)-VX(v2))^2 + (VY(v1)-VY(v2))^2 );

                    % find find volume node numbers of left and right nodes 
                    vidM = vmapM(flist1, k1); vidP = vmapM(flist2, k2);
                    x1 = obj.x(vidM); y1 = obj.y(vidM); x2 = obj.x(vidP); y2 = obj.y(vidP);
                    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

                    % Compute distance matrix
                    D = (x1 -x2').^2 + (y1-y2').^2;
                    [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
                    vmapP(flist1(idM), k1) = vidP(idP);
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