classdef tri_mesh < ndg_lib.mesh.mesh2d
    %TRI_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = tri_mesh(cell, varargin)
            obj = obj@ndg_lib.mesh.mesh2d(cell, varargin{:});
            
            % check input element
            if( ne(obj.cell.type, ndg_lib.std_cell_type.Tri) )
                error(['Input cell type ', cell.type, 'is not triangle!'])
            end
            
        end% func
        
        obj = refine(obj, multi_rate); % ¼ÓÃÜÍø¸ñ
    end% methods
    
end

