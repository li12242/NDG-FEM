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
            if ~isempty(BC)
                obj.mapO = MultiRegions.getBoundMapList(obj, BC, OUTFLOW);
                obj.mapI = MultiRegions.getBoundMapList(obj, BC, INFLOW);
                obj.mapW = MultiRegions.getBoundMapList(obj, BC, WALL);
            end% if
        end% function
    end% methods
end% classedf