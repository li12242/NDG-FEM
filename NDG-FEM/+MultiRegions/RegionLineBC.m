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
        obj.mapO = MultiRegions.GetBoundMapList(obj,BC,OUTFLOW);
        obj.mapI = MultiRegions.GetBoundMapList(obj,BC,INFLOW);
        obj.mapW = MultiRegions.GetBoundMapList(obj,BC,WALL);
    end
end% methods public
end %classdef Region1D
