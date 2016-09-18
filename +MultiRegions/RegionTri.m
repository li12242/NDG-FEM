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
            
            [obj.EToE, obj.EToF,  obj] = MultiRegions.Connect2D(obj,EToV);
            [obj.vmapM,obj.vmapP, obj] = MultiRegions.BuildMap...
                (obj, VX, VY, EToV, obj.EToE, obj.EToF);
            
            obj.fScale = obj.sJ./obj.J(obj.vmapM);
        end% function
    end% methods
end% classedf