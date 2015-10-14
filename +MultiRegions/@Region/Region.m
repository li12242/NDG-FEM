classdef Region < handle
    % physical meshes
    % Region properties:
    %   Shape    - base element type
    %   nElement - No. of element
    %   nNode    - No. of node
    %   vmapM    - face list to node list (in mesh)
    %   vmapP    - face list to neighbour node list (in mesh)
    %   sJ       - face factor: surface element jacobi factor
    %   J        - element factor: element Jacobi factor
    %   EToE     - element to element
    %   EToF     - element to face
    % Region properties(private):
    %   SpFToV   - face to vertices
    % Region methods: none
    % 
    properties(SetAccess = protected)
        Shape       % base element type
        nElement    % No. of element
        nNode       % No. of node
        vmapM       % face list to node list (in mesh)
        vmapP       % face list to neighbour node list (in mesh)
        sJ      % face factor: surface Jacobi
        J       % element factor: Jacobi transfer factor
        EToE    % element to element
        EToF    % element to face
    end
    properties(GetAccess=protected)
        SpFToV
    end% properties
    
    methods
        function obj = Region(shape, EToV)
            % constructor
            obj.Shape = shape;
            obj.nElement = size(EToV, 1);
            obj.nNode = obj.nElement*shape.nNode;
        end
    end
end