classdef quad_mesh < ndg_lib.mesh.mesh2d
    %QUAD_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = quad_mesh(cell, varargin)
            obj = obj@ndg_lib.mesh.mesh2d(cell, varargin{:});
            
            % check input element
            if( ne(obj.cell.type, ndg_lib.std_cell_type.Quad) )
                error(['Input cell type ', cell.type, 'is not quadrilateral!'])
            end
            
        end% func
    end% methods
    
end

