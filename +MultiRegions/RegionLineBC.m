classdef RegionLineBC < MultiRegions.RegionLine
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
%   x        - node coordinate
%   rx       - element factor: dr/dx
%   nx       - face factor: normal vector
%   fScale   - face factor: sJ/J(surface node)
% RegionLine properties(private):
%   SpFToV   - face to vertices
% RegionLineBC properties(public):
%   mapI     - Inflow  BC node to face list 
%   mapO     - Outflow BC node to face list
% RegionLineBC methods: none
% Usages:
%   mesh = RegionLineBC(line, EToV, VX, BC)
% Input:    
%   VX  - vertice coordinate, size [nVertic  x1]
%   EToV- element to vertice, size [nElement x2]
%   BC  - boundary types, [BCtypeID, verticeID]
properties       
    mapI    % Inflow  BC node to face list 
    mapO    % Outflow BC node to face list
end

methods
    %% RegionLineBC
    % Construction function of line mesh with boundary points.
    % Input:    
    %   EToV - element to vertice, size [nElementx2]
    %   VX - vertice coordinate, size [nVerticx1]
    %   VY - vertice coordinate, size [nVerticx1]
    %   BC - boundary types, size [BCtypeID, nodeId]
    % Usages:
    %
    function obj = RegionLineBC(line, EToV, VX, BC)
        obj = obj@MultiRegions.RegionLine(line, EToV, VX);

        MultiRegions.BCType(); % include boundary type ID
        obj.mapO = getBoundMapList(obj,BC,OUTFLOW);
        obj.mapI = getBoundMapList(obj,BC,INFLOW);
        obj.mapW = getBoundMapList(obj,BC,WALL);
    end
end% methods public

methods(Hidden)
    %% getBoundMapList
    % 
    %
    function map = getBoundMapList(obj,BC,phyID)
        bclist = (BC(:,1)==phyID);
        partialBC = BC(bclist, :);
        nBC = size(partialBC,1); 
        Nv = size(obj.SpFToV, 2); % number of vertice
        VToBF = spalloc(Nv, nBC, 2*nBC);
        for ib = 1:nBC
            VToBF(partialBC(ib,2),ib) = 1;
        end
        FToBF = obj.SpFToV*VToBF;
        [faces1, ~] = find(FToBF ==1);

        Nfaces = 2;
        element1 = floor( (faces1-1)/Nfaces )  + 1; 
        face1    =   mod( (faces1-1), Nfaces ) + 1;

        Nfp = obj.Shape.nFaceNode/obj.Shape.nFace; 
        map = zeros(Nfp, nBC);
        for i=1:nBC
            localfacelist = obj.Shape.getFaceListAtFace(face1(i));
            facelist = (element1(i)-1)*obj.Shape.nFaceNode ...
                + localfacelist;
            map(:,i) = facelist;
        end
    end% function

end %methods private
end %classdef Region1D
