classdef RegionQuadBC < MultiRegions.RegionQuad
    % physical meshes with Quadrilateral element
    % RegionQuad properties(inherit):
    %   Shape    - base element type
    %   nElement - No. of element
    %   nNode    - No. of node
    %   vmapM    - face list to node list (in mesh)
    %   vmapP    - face list to neighbour node list (in mesh)
    %   sJ       - face factor: surface element Jacobi
    %   J        - element factor: Jacobi transfer factor
    %   EToE     - element to element
    %   EToF     - element to face
    % RegionQuad properties:
    %   x       - node coordinate
    %   y       - node coordinate
    %   rx      - element factor: dr/dx
    %   sx      - element factor: ds/dx
    %   ry      - element factor: dr/dy
    %   sy      - element factor: ds/dy
    %   nx      - face factor: normal vector
    %   ny      - face factor: normal vector
    %   fScale  - face factor: sJ/J(surface node)
    % RegionQuad methods: none
    properties
        mapO    % bc(outflow) list to face list
        mapI    % bc(inflow) list to face list
        mapW    % bc(wall) list to face list
    end% properties
    methods
        function obj = RegionQuadBC(shape, EToV, VX, VY, BC)
            % Input:    EToV - element to vertice, size [nElementx3]
            %           VX - vertice coordinate, size [nVerticx1]
            %           VY - vertice coordinate, size [nVerticx1]
            %           BC - boundary types, [BCtypeID, node1, node2]
            obj = obj@MultiRegions.RegionQuad(shape, EToV, VX, VY);
            
            MultiRegions.BCType
            if ~isempty(BC)
                obj.mapO = MultiRegions.GetBoundMapList(obj, BC, OUTFLOW);
                obj.mapI = MultiRegions.GetBoundMapList(obj, BC, INFLOW);
                obj.mapW = MultiRegions.GetBoundMapList(obj, BC, WALL);
            end% if
        end% function
    end% methods
end% class