function [ limiter ] = BJ_wrapper( mesh )
%BJ_WRAPPER Construct the BJ slope limiter for different types of meshes.
%   Detailed explanation goes here

switch mesh.cell.type
    case ndg_lib.std_cell_type.Line
        limiter = ndg_utility.limiter.BJ.BJ_line(mesh);
    case ndg_lib.std_cell_type.Tri
    case ndg_lib.std_cell_type.Quad
end

end
