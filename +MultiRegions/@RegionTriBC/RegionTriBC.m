classdef RegionTriBC < MultiRegions.RegionTri
    % physical meshes with Line element
    % RegionTriBC properties(inherit):
    %   Shape    - base element type
    %   nElement - No. of element
    %   nNode    - No. of node
    %   vmapM    - face list to node list (in mesh)
    %   vmapP    - face list to neighbour node list (in mesh)
    %   sJ       - face factor: surface element Jacobi
    %   J        - element factor: Jacobi transfer factor
    %   EToE     - element to element
    %   EToF     - element to face
    %   x       - node coordinate
    %   y       - node coordinate
    %   rx      - element factor: dr/dx
    %   sx      - element factor: ds/dx
    %   ry      - element factor: dr/dy
    %   sy      - element factor: ds/dy
    %   nx      - face factor: normal vector
    %   ny      - face factor: normal vector
    %   fScale  - face factor: sJ/J(surface node)
    % RegionTriBC properties(inherit, private):
    %   SpFToV  - face to vertices
    % RegionTriBC properties:
    %   mapO    - bc(outflow) list to face list, size [Nfp x nBC]
    %   mapW    - bc(wall) list to face list [Nfp x nBC]
    % RegionTri methods: none
    properties
        mapO    % bc(outflow) list to face list
        mapI    % bc(inflow) list to face list
        mapW    % bc(wall) list to face list
    end% properties
    methods
        function obj = RegionTriBC(shape, EToV, VX, VY, BC)
            % Input:    EToV - element to vertice, size [nElementx3]
            %           VX - vertice coordinate, size [nVerticx1]
            %           VY - vertice coordinate, size [nVerticx1]
            %           BC - boundary types, [BCtypeID, node1, node2]
            obj = obj@MultiRegions.RegionTri(shape, EToV, VX, VY);
            
            MultiRegions.BCType
            obj.mapO = getBoundMapList(obj, BC, OUTFLOW);
            obj.mapI = getBoundMapList(obj, BC, INFLOW);
            obj.mapW = getBoundMapList(obj, BC, WALL);
        end% function
    end% methods
    
    methods(Hidden)
        function map = getBoundMapList(obj,BC,phyID)
            % Input:    phyID - physical ID
            bclist = find(BC(:,1)==phyID);
            partialBC = BC(bclist, :);
            nBC = size(partialBC,1); Nv = size(obj.SpFToV, 2);
            VToBF = spalloc(Nv, nBC, 2*nBC);
            for ib = 1:nBC
                VToBF(partialBC(ib,[2,3]),ib) = 1;
            end
            FToBF = obj.SpFToV*VToBF;
            [faces1, ~] = find(FToBF == 2);
            
            Nfaces = 3;
            element1 = floor( (faces1-1)/Nfaces )  + 1; 
            face1    =   mod( (faces1-1), Nfaces ) + 1;
            
            Nfp = obj.Shape.nFaceNode/obj.Shape.nFace; 
            map = zeros(Nfp, nBC);
            for i=1:nBC
                localfacelist = obj.Shape.getFaceListAtFace(face1(i));
                facelist = (element1(i)-1)*obj.Shape.nFaceNode + localfacelist;
                map(:,i) = facelist;
            end
        end% function
        
    end% methods
end% classedf